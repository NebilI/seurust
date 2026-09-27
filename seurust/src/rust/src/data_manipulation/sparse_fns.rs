use crate::sparse::{
    build_dgcmatrix, rmatrix_from_column_major, CscSlots, CscView, CsrView, RowIndex,
};
use extendr_api::prelude::*;
use extendr_api::GetSexp;
use extendr_ffi::Rf_runif;
use rayon::prelude::*;
use std::collections::HashMap;

fn log1p(x: f64) -> f64 {
    x.ln_1p()
}

fn expm1(x: f64) -> f64 {
    x.exp_m1()
}

/// Normalize CSC `x` values in place.
pub fn log_norm_impl(
    x: &mut [f64],
    p: &[i32],
    col_sums: &[f64],
    ncols: usize,
    scale_factor: i32,
    _display_progress: bool,
) {
    let scale = scale_factor as f64;

    for col in 0..ncols {
        for idx in p[col] as usize..p[col + 1] as usize {
            x[idx] = log1p(x[idx] / col_sums[col] * scale);
        }
    }
}

pub fn log_norm_owned_impl(mat: &mut CscSlots, scale_factor: i32, display_progress: bool) {
    let view = CscView {
        x: &mat.x,
        i: &mat.i,
        p: &mat.p,
        nrows: mat.nrows,
        ncols: mat.ncols,
    };
    let col_sums = view.col_sums();
    let ncols = mat.ncols as usize;
    log_norm_impl(
        &mut mat.x,
        &mat.p,
        &col_sums,
        ncols,
        scale_factor,
        display_progress,
    );
}

pub fn run_umi_sampling_impl(mut mat: CscSlots, sample_val: i32, upsample: bool) -> CscSlots {
    let col_sums = mat.col_sums();
    let target = sample_val as f64;
    let ncols = mat.ncols as usize;

    for col in 0..ncols {
        let col_sum = col_sums[col];
        if upsample || col_sum > target {
            for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
                let mut entry = mat.x[idx] * target / col_sum;
                let frac = entry.fract();
                if frac != 0.0 {
                    let rn = unsafe { Rf_runif(0.0, 1.0) };
                    entry = if frac <= rn {
                        entry.floor()
                    } else {
                        entry.ceil()
                    };
                }
                mat.x[idx] = entry;
            }
        }
    }

    mat
}

pub fn run_umi_sampling_per_cell_impl(
    mut mat: CscSlots,
    sample_val: &[f64],
    upsample: bool,
) -> CscSlots {
    let col_sums = mat.col_sums();
    let ncols = mat.ncols as usize;

    for col in 0..ncols {
        let col_sum = col_sums[col];
        let target = sample_val[col];
        if upsample || col_sum > target {
            for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
                let mut entry = mat.x[idx] * target / col_sum;
                let frac = entry.fract();
                if frac != 0.0 {
                    let rn = unsafe { Rf_runif(0.0, 1.0) };
                    entry = if frac <= rn {
                        entry.floor()
                    } else {
                        entry.ceil()
                    };
                }
                mat.x[idx] = entry;
            }
        }
    }

    mat
}

/// R interns CHARSXP values, so equal names share a pointer.
fn char_key(name: &Rstr) -> usize {
    unsafe { name.get() as usize }
}

pub fn row_merge_matrices_impl(
    mat1: CsrView<'_>,
    mat2: CsrView<'_>,
    mat1_rownames: &Strings,
    mat2_rownames: &Strings,
    all_rownames: &Strings,
) -> CscSlots {
    let mut mat1_map: HashMap<usize, usize> = HashMap::with_capacity(mat1_rownames.len());
    for (idx, name) in mat1_rownames.iter().enumerate() {
        mat1_map.insert(char_key(name), idx);
    }
    let mut mat2_map: HashMap<usize, usize> = HashMap::with_capacity(mat2_rownames.len());
    for (idx, name) in mat2_rownames.iter().enumerate() {
        mat2_map.insert(char_key(name), idx);
    }

    let num_rows = all_rownames.len();
    let num_col1 = mat1.ncols as usize;
    let num_cols = num_col1 + mat2.ncols as usize;

    // usize::MAX marks a name missing from that input.
    let mut src1 = vec![usize::MAX; num_rows];
    let mut src2 = vec![usize::MAX; num_rows];
    for (row, name) in all_rownames.iter().enumerate() {
        let key = char_key(name);
        if let Some(&src) = mat1_map.get(&key) {
            src1[row] = src;
        }
        if let Some(&src) = mat2_map.get(&key) {
            src2[row] = src;
        }
    }

    // Count nonzeros per output column, then scatter. CSR rows are addressed
    // from `p` in constant time, and rows are filled in order so each column
    // ends up sorted without a triplet sort.
    let mut col_counts = vec![0usize; num_cols];
    let mut nnz = 0usize;
    for row in 0..num_rows {
        if src1[row] != usize::MAX {
            let start = mat1.p[src1[row]] as usize;
            let end = mat1.p[src1[row] + 1] as usize;
            nnz += end - start;
            for idx in start..end {
                col_counts[mat1.j[idx] as usize] += 1;
            }
        }
        if src2[row] != usize::MAX {
            let start = mat2.p[src2[row]] as usize;
            let end = mat2.p[src2[row] + 1] as usize;
            nnz += end - start;
            for idx in start..end {
                col_counts[num_col1 + mat2.j[idx] as usize] += 1;
            }
        }
    }

    let mut p = vec![0i32; num_cols + 1];
    let mut cursor = vec![0usize; num_cols];
    let mut acc = 0usize;
    for col in 0..num_cols {
        p[col] = acc as i32;
        cursor[col] = acc;
        acc += col_counts[col];
    }
    p[num_cols] = acc as i32;

    let mut i = vec![0i32; nnz];
    let mut x = vec![0.0; nnz];
    for row in 0..num_rows {
        let row_i = row as i32;
        if src1[row] != usize::MAX {
            let start = mat1.p[src1[row]] as usize;
            let end = mat1.p[src1[row] + 1] as usize;
            for idx in start..end {
                let col = mat1.j[idx] as usize;
                let dest = cursor[col];
                cursor[col] = dest + 1;
                i[dest] = row_i;
                x[dest] = mat1.x[idx];
            }
        }
        if src2[row] != usize::MAX {
            let start = mat2.p[src2[row]] as usize;
            let end = mat2.p[src2[row] + 1] as usize;
            for idx in start..end {
                let col = num_col1 + mat2.j[idx] as usize;
                let dest = cursor[col];
                cursor[col] = dest + 1;
                i[dest] = row_i;
                x[dest] = mat2.x[idx];
            }
        }
    }

    CscSlots {
        x,
        i,
        p,
        nrows: num_rows as i32,
        ncols: num_cols as i32,
    }
}

fn scale_clip(value: f64, scale_max: f64) -> f64 {
    if value > scale_max {
        scale_max
    } else {
        value
    }
}

fn gene_mean_sdev(
    view: &CscView<'_>,
    row_index: &RowIndex,
    gene: usize,
    n_cells: usize,
    scale: bool,
    center: bool,
) -> (f64, f64) {
    let range = row_index.row_range(gene);
    let col_mean: f64 = (range.start..range.end)
        .map(|pos| view.x[row_index.row_x_idx[pos]])
        .sum::<f64>()
        / n_cells as f64;

    let mut col_sdev = 1.0;
    if scale {
        let nn_zero = range.len();
        let mut var_sum = 0.0;
        if center {
            for pos in range.clone() {
                let val = view.x[row_index.row_x_idx[pos]];
                var_sum += (val - col_mean).powi(2);
            }
            var_sum += col_mean.powi(2) * (n_cells - nn_zero) as f64;
        } else {
            var_sum = (range.start..range.end)
                .map(|pos| view.x[row_index.row_x_idx[pos]].powi(2))
                .sum();
        }
        col_sdev = (var_sum / (n_cells - 1) as f64).sqrt();
    }

    let mean = if center { col_mean } else { 0.0 };
    (mean, col_sdev)
}

pub fn fast_sparse_row_scale_impl(
    view: CscView<'_>,
    scale: bool,
    center: bool,
    scale_max: f64,
    _display_progress: bool,
) -> RMatrix<f64> {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    let row_index = RowIndex::from_csc_view(&view);
    let mut data = vec![0.0; n_genes * n_cells];
    let out_addr = data.as_mut_ptr() as usize;

    (0..n_genes).into_par_iter().for_each(|gene| {
        let (mean, col_sdev) = gene_mean_sdev(&view, &row_index, gene, n_cells, scale, center);
        let inv_sdev = 1.0 / col_sdev;
        let out_ptr = out_addr as *mut f64;
        unsafe {
            for cell in 0..n_cells {
                let value = scale_clip((0.0 - mean) * inv_sdev, scale_max);
                *out_ptr.add(gene + cell * n_genes) = value;
            }
            for pos in row_index.row_range(gene) {
                let cell = row_index.row_cols[pos];
                let val = view.x[row_index.row_x_idx[pos]];
                let value = scale_clip((val - mean) * inv_sdev, scale_max);
                *out_ptr.add(gene + cell * n_genes) = value;
            }
        }
    });

    rmatrix_from_column_major(&data, n_genes, n_cells)
}

pub fn fast_sparse_row_scale_with_known_stats_impl(
    view: CscView<'_>,
    mu: &[f64],
    sigma: &[f64],
    scale: bool,
    center: bool,
    scale_max: f64,
    _display_progress: bool,
) -> RMatrix<f64> {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    let row_index = RowIndex::from_csc_view(&view);
    let mut data = vec![0.0; n_genes * n_cells];
    let out_addr = data.as_mut_ptr() as usize;

    (0..n_genes).into_par_iter().for_each(|gene| {
        let col_mean = if center { mu[gene] } else { 0.0 };
        let col_sdev = if scale { sigma[gene] } else { 1.0 };
        let inv_sdev = 1.0 / col_sdev;
        let out_ptr = out_addr as *mut f64;
        unsafe {
            for cell in 0..n_cells {
                let value = scale_clip((0.0 - col_mean) * inv_sdev, scale_max);
                *out_ptr.add(gene + cell * n_genes) = value;
            }
            for pos in row_index.row_range(gene) {
                let cell = row_index.row_cols[pos];
                let val = view.x[row_index.row_x_idx[pos]];
                let value = scale_clip((val - col_mean) * inv_sdev, scale_max);
                *out_ptr.add(gene + cell * n_genes) = value;
            }
        }
    });

    rmatrix_from_column_major(&data, n_genes, n_cells)
}

/// Parallel per-row sums over CSC columns.
///
/// Nonzeros of one gene are not contiguous in a `dgCMatrix`, but each entry
/// can still be added into a dense row accumulator while the column data is
/// scanned in order. That avoids building a row index (and the unused column-id
/// buffer that went with it) and keeps the hot loop on sequential `x`/`i`.
fn par_row_sums(view: &CscView<'_>, step: &(impl Fn(usize, f64, &mut f64) + Sync)) -> Vec<f64> {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    (0..n_cells)
        .into_par_iter()
        .fold(
            || vec![0.0; n_genes],
            |mut acc, col| {
                for idx in view.p[col] as usize..view.p[col + 1] as usize {
                    let row = view.i[idx] as usize;
                    step(row, view.x[idx], &mut acc[row]);
                }
                acc
            },
        )
        .reduce(
            || vec![0.0; n_genes],
            |mut left, right| {
                for (dst, src) in left.iter_mut().zip(right) {
                    *dst += src;
                }
                left
            },
        )
}

/// Like [`par_row_sums`], also counting how many stored entries each row has.
fn par_row_sums_and_counts(
    view: &CscView<'_>,
    step: &(impl Fn(usize, f64, &mut f64, &mut usize) + Sync),
) -> (Vec<f64>, Vec<usize>) {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    (0..n_cells)
        .into_par_iter()
        .fold(
            || (vec![0.0; n_genes], vec![0usize; n_genes]),
            |mut acc, col| {
                for idx in view.p[col] as usize..view.p[col + 1] as usize {
                    let row = view.i[idx] as usize;
                    step(row, view.x[idx], &mut acc.0[row], &mut acc.1[row]);
                }
                acc
            },
        )
        .reduce(
            || (vec![0.0; n_genes], vec![0usize; n_genes]),
            |mut left, right| {
                for gene in 0..n_genes {
                    left.0[gene] += right.0[gene];
                    left.1[gene] += right.1[gene];
                }
                left
            },
        )
}

fn fast_exp_mean_values(view: &CscView<'_>) -> Vec<f64> {
    let ncols_f = view.ncols as f64;
    par_row_sums(view, &|_, val, sum| {
        *sum += expm1(val);
    })
    .into_iter()
    .map(|sum| log1p(sum / ncols_f))
    .collect()
}

pub fn fast_exp_mean_impl(view: CscView<'_>, _display_progress: bool) -> Doubles {
    Doubles::from_values(fast_exp_mean_values(&view))
}

pub fn sparse_row_var2_impl(view: CscView<'_>, mu: &[f64], _display_progress: bool) -> Doubles {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    let denom = (n_cells as f64) - 1.0;

    let (sums, nnz_rows): (Vec<f64>, Vec<usize>) = (0..n_cells)
        .into_par_iter()
        .fold(
            || (vec![0.0; n_genes], vec![0usize; n_genes]),
            |mut acc, col| {
                for idx in view.p[col] as usize..view.p[col + 1] as usize {
                    let row = view.i[idx] as usize;
                    acc.1[row] += 1;
                    let diff = view.x[idx] - mu[row];
                    acc.0[row] += diff * diff;
                }
                acc
            },
        )
        .reduce(
            || (vec![0.0; n_genes], vec![0usize; n_genes]),
            |mut left, right| {
                for gene in 0..n_genes {
                    left.0[gene] += right.0[gene];
                    left.1[gene] += right.1[gene];
                }
                left
            },
        );

    let all_vars: Vec<f64> = (0..n_genes)
        .map(|gene_idx| {
            let n_zero = n_cells - nnz_rows[gene_idx];
            let mu_i = mu[gene_idx];
            (sums[gene_idx] + mu_i * mu_i * n_zero as f64) / denom
        })
        .collect();

    Doubles::from_values(all_vars)
}

fn sparse_row_var_std_values(view: &CscView<'_>, mu: &[f64], sd: &[f64], vmax: f64) -> Vec<f64> {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    let denom = n_cells as f64 - 1.0;
    // sd == 0 matches the C++ `continue`: the row stays 0 and is not standardized.
    let (sums, nnz_rows) = par_row_sums_and_counts(view, &|row, val, sum, count| {
        if sd[row] == 0.0 {
            return;
        }
        *count += 1;
        let standardized = ((val - mu[row]) / sd[row]).min(vmax);
        *sum += standardized * standardized;
    });

    (0..n_genes)
        .map(|gene| {
            if sd[gene] == 0.0 {
                return 0.0;
            }
            let n_zero = n_cells - nnz_rows[gene];
            let zero = (0.0 - mu[gene]) / sd[gene];
            (sums[gene] + zero * zero * n_zero as f64) / denom
        })
        .collect()
}

pub fn sparse_row_var_std_impl(
    view: CscView<'_>,
    mu: &[f64],
    sd: &[f64],
    vmax: f64,
    _display_progress: bool,
) -> Doubles {
    Doubles::from_values(sparse_row_var_std_values(&view, mu, sd, vmax))
}

fn fast_log_vmr_values(view: &CscView<'_>) -> Vec<f64> {
    let nrows = view.nrows as usize;
    let ncols = view.ncols as usize;
    let ncols_f = ncols as f64;
    let (sums, nnz) = par_row_sums_and_counts(view, &|_, val, sum, count| {
        *count += 1;
        *sum += expm1(val);
    });
    let rm: Vec<f64> = sums.into_iter().map(|sum| sum / ncols_f).collect();
    let sq = par_row_sums(view, &|row, val, sum| {
        let diff = expm1(val) - rm[row];
        *sum += diff * diff;
    });

    (0..nrows)
        .map(|row| {
            let v = (sq[row] + (ncols - nnz[row]) as f64 * rm[row] * rm[row]) / (ncols_f - 1.0);
            (v / rm[row]).ln()
        })
        .collect()
}

pub fn fast_log_vmr_impl(view: CscView<'_>, _display_progress: bool) -> Doubles {
    Doubles::from_values(fast_log_vmr_values(&view))
}

fn sparse_row_var_values(view: &CscView<'_>) -> Vec<f64> {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    let n_cells_f = n_cells as f64;
    let denom = n_cells_f - 1.0;
    let (sums, nnz) = par_row_sums_and_counts(view, &|_, val, sum, count| {
        *count += 1;
        *sum += val;
    });
    let rm: Vec<f64> = sums.into_iter().map(|sum| sum / n_cells_f).collect();
    let sq = par_row_sums(view, &|row, val, sum| {
        let diff = val - rm[row];
        *sum += diff * diff;
    });

    (0..n_genes)
        .map(|row| (sq[row] + (n_cells - nnz[row]) as f64 * rm[row] * rm[row]) / denom)
        .collect()
}

pub fn sparse_row_var_impl(view: CscView<'_>, _display_progress: bool) -> Doubles {
    Doubles::from_values(sparse_row_var_values(&view))
}

pub fn replace_cols_impl(
    mat: CscView<'_>,
    col_idx: &[i32],
    replacement: CscView<'_>,
) -> extendr_api::Result<Robj> {
    let ncols = mat.ncols as usize;
    // Last duplicate index wins, matching sequential column assignment.
    let mut which = vec![-1i32; ncols];
    for (rep_idx, &col) in col_idx.iter().enumerate() {
        if col >= 0 {
            let col = col as usize;
            if col < ncols {
                which[col] = rep_idx as i32;
            }
        }
    }

    let mut nnz = 0usize;
    for col in 0..ncols {
        let (p, src_col) = if which[col] >= 0 {
            (replacement.p, which[col] as usize)
        } else {
            (mat.p, col)
        };
        nnz += (p[src_col + 1] - p[src_col]) as usize;
    }

    build_dgcmatrix(mat.nrows, mat.ncols, nnz, |x, i_out, p_out| {
        let mut dest = 0usize;
        for col in 0..ncols {
            p_out[col] = dest as i32;
            let (src_i, src_x, start, end) = if which[col] >= 0 {
                let rep = which[col] as usize;
                let start = replacement.p[rep] as usize;
                let end = replacement.p[rep + 1] as usize;
                (replacement.i, replacement.x, start, end)
            } else {
                let start = mat.p[col] as usize;
                let end = mat.p[col + 1] as usize;
                (mat.i, mat.x, start, end)
            };
            let n = end - start;
            i_out[dest..dest + n].copy_from_slice(&src_i[start..end]);
            x[dest..dest + n].copy_from_slice(&src_x[start..end]);
            dest += n;
        }
        p_out[ncols] = dest as i32;
    })
}

/// Per-row neighbor lists of a graph (cells x cells), ordered by increasing
/// distance. Row `k` of the outputs holds the 1-based column indices and
/// values of row `k` of `mat`, as in Seurat's `GraphToNeighborHelper`.
pub fn graph_to_neighbor_lists(mat: &CscSlots) -> Result<(usize, Vec<f64>, Vec<f64>), String> {
    let csr = mat.to_cs_mat().to_csr();
    let nrows = csr.rows();
    let n_neighbors = csr
        .outer_iterator()
        .next()
        .map(|row| row.nnz())
        .unwrap_or(0);

    // Column-major nrows x n_neighbors buffers, ready to hand to R.
    let mut nn_idx = vec![0.0; nrows * n_neighbors];
    let mut nn_dist = vec![0.0; nrows * n_neighbors];

    for (k, row) in csr.outer_iterator().enumerate() {
        if row.nnz() != n_neighbors {
            return Err(format!(
                "Not all cells have an equal number of neighbors \
                 (row 1 has {n_neighbors}, row {} has {}).",
                k + 1,
                row.nnz()
            ));
        }

        let row_idx = row.indices();
        let row_dist = row.data();
        let mut order: Vec<usize> = (0..row_dist.len()).collect();
        order.sort_by(|&a, &b| {
            row_dist[a]
                .partial_cmp(&row_dist[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        for (i, &ord) in order.iter().enumerate() {
            nn_idx[k + i * nrows] = (row_idx[ord] + 1) as f64;
            nn_dist[k + i * nrows] = row_dist[ord];
        }
    }

    Ok((n_neighbors, nn_idx, nn_dist))
}

pub fn graph_to_neighbor_helper_impl(mat: CscSlots) -> extendr_api::Result<Robj> {
    let nrows = mat.nrows as usize;
    let (n_neighbors, nn_idx, nn_dist) =
        graph_to_neighbor_lists(&mat).map_err(extendr_api::Error::Other)?;

    Ok(Robj::from(List::from_values([
        Robj::from(rmatrix_from_column_major(&nn_idx, nrows, n_neighbors)),
        Robj::from(rmatrix_from_column_major(&nn_dist, nrows, n_neighbors)),
    ])))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::CscSlots;

    fn toy_csc() -> CscSlots {
        CscSlots {
            x: vec![1.0, 2.0, 3.0],
            i: vec![0, 2, 1],
            p: vec![0, 1, 2, 3],
            nrows: 3,
            ncols: 3,
        }
    }

    fn toy_view<'a>(mat: &'a CscSlots) -> CscView<'a> {
        CscView {
            x: &mat.x,
            i: &mat.i,
            p: &mat.p,
            nrows: mat.nrows,
            ncols: mat.ncols,
        }
    }

    #[test]
    fn log_norm_scales_columns() {
        let mut mat = toy_csc();
        log_norm_owned_impl(&mut mat, 10_000, false);
        assert!(mat.x.iter().all(|v| v.is_finite() && *v >= 0.0));
    }

    #[test]
    fn graph_to_neighbor_lists_reads_rows_sorted_by_distance() {
        // Asymmetric 3x3 graph, dense row-major:
        //   row 0: [0.0, 0.5, 0.2]
        //   row 1: [0.9, 0.0, 0.1]
        //   row 2: [0.3, 0.4, 0.0]
        let mat = CscSlots {
            x: vec![0.9, 0.3, 0.5, 0.4, 0.2, 0.1],
            i: vec![1, 2, 0, 2, 0, 1],
            p: vec![0, 2, 4, 6],
            nrows: 3,
            ncols: 3,
        };
        let (k, idx, dist) = graph_to_neighbor_lists(&mat).unwrap();
        assert_eq!(k, 2);
        // Column-major 3 x 2 outputs.
        assert_eq!(idx, vec![3.0, 3.0, 1.0, 2.0, 1.0, 2.0]);
        assert_eq!(dist, vec![0.2, 0.1, 0.3, 0.5, 0.9, 0.4]);
    }

    #[test]
    fn graph_to_neighbor_lists_rejects_ragged_rows() {
        let mat = CscSlots {
            x: vec![1.0, 1.0, 1.0],
            i: vec![1, 0, 1],
            p: vec![0, 1, 3],
            nrows: 2,
            ncols: 2,
        };
        assert!(graph_to_neighbor_lists(&mat).is_err());
    }

    /// Sequential column scan of the C++ formulas, used as the reference.
    fn reference_stats(mat: &CscSlots) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let n_genes = mat.nrows as usize;
        let n_cells = mat.ncols as usize;
        let n_cells_f = n_cells as f64;
        let mut exp_sum = vec![0.0; n_genes];
        let mut raw_sum = vec![0.0; n_genes];
        let mut nnz = vec![0usize; n_genes];
        for col in 0..n_cells {
            for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
                let row = mat.i[idx] as usize;
                let val = mat.x[idx];
                exp_sum[row] += val.exp_m1();
                raw_sum[row] += val;
                nnz[row] += 1;
            }
        }
        let exp_mean: Vec<f64> = exp_sum
            .iter()
            .map(|sum| (sum / n_cells_f).ln_1p())
            .collect();
        let rm: Vec<f64> = raw_sum.iter().map(|sum| sum / n_cells_f).collect();
        let exp_rm: Vec<f64> = exp_sum.iter().map(|sum| sum / n_cells_f).collect();
        let mut var = vec![0.0; n_genes];
        let mut vmr = vec![0.0; n_genes];
        let mut std_var = vec![0.0; n_genes];
        let mu = [0.5, 1.0];
        let sd = [2.0, 0.0];
        let vmax = 0.5;
        for col in 0..n_cells {
            for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
                let row = mat.i[idx] as usize;
                let val = mat.x[idx];
                let d = val - rm[row];
                var[row] += d * d;
                let e = val.exp_m1() - exp_rm[row];
                vmr[row] += e * e;
                if sd[row] != 0.0 {
                    let standardized = ((val - mu[row]) / sd[row]).min(vmax);
                    std_var[row] += standardized * standardized;
                }
            }
        }
        for row in 0..n_genes {
            let n_zero = n_cells - nnz[row];
            var[row] = (var[row] + n_zero as f64 * rm[row] * rm[row]) / (n_cells_f - 1.0);
            let v = (vmr[row] + n_zero as f64 * exp_rm[row] * exp_rm[row]) / (n_cells_f - 1.0);
            vmr[row] = (v / exp_rm[row]).ln();
            if sd[row] == 0.0 {
                std_var[row] = 0.0;
            } else {
                let zero = (0.0 - mu[row]) / sd[row];
                std_var[row] = (std_var[row] + zero * zero * n_zero as f64) / (n_cells_f - 1.0);
            }
        }
        (exp_mean, var, vmr, std_var)
    }

    #[test]
    fn row_stats_match_column_scan_formula() {
        // row 0: 1, 0, 3, 0    row 1: 0, 4, 0, 0
        let mat = CscSlots {
            x: vec![1.0, 4.0, 3.0],
            i: vec![0, 1, 0],
            p: vec![0, 1, 2, 3, 3],
            nrows: 2,
            ncols: 4,
        };
        let view = toy_view(&mat);
        let (exp_mean, var, vmr, std_var) = reference_stats(&mat);
        let got_exp = fast_exp_mean_values(&view);
        let got_var = sparse_row_var_values(&view);
        let got_vmr = fast_log_vmr_values(&view);
        let got_std = sparse_row_var_std_values(&view, &[0.5, 1.0], &[2.0, 0.0], 0.5);
        for (got, expected) in [
            (&got_exp, &exp_mean),
            (&got_var, &var),
            (&got_vmr, &vmr),
            (&got_std, &std_var),
        ] {
            assert_eq!(got.len(), expected.len());
            for (g, e) in got.iter().zip(expected) {
                if e.is_nan() {
                    assert!(g.is_nan());
                } else {
                    assert!((g - e).abs() <= 1e-12, "{g} vs {e}");
                }
            }
        }
        // sd == 0 leaves that row at 0, including when the row has nonzeros.
        assert_eq!(got_std[1], 0.0);
    }

    #[test]
    fn fast_sparse_row_scale_computes_finite_stats() {
        let mat = toy_csc();
        let view = toy_view(&mat);
        let row_index = RowIndex::from_csc_view(&view);
        let (mean, sdev) = gene_mean_sdev(&view, &row_index, 0, 3, true, true);
        assert!(mean.is_finite() && sdev.is_finite() && sdev > 0.0);
    }
}
