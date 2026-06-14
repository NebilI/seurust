use crate::sparse::{
    csc_from_triplets, rmatrix_from_column_major, CscSlots, CscView, CsrSlots, RowIndex,
};
use extendr_api::prelude::*;
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

pub fn row_merge_matrices_impl(
    mat1: CsrSlots,
    mat2: CsrSlots,
    mat1_rownames: &[String],
    mat2_rownames: &[String],
    all_rownames: &[String],
) -> CscSlots {
    let mat1 = mat1.to_cs_mat();
    let mat2 = mat2.to_cs_mat();

    let mut mat1_map: HashMap<&str, usize> = HashMap::new();
    for (idx, name) in mat1_rownames.iter().enumerate() {
        mat1_map.insert(name.as_str(), idx);
    }
    let mut mat2_map: HashMap<&str, usize> = HashMap::new();
    for (idx, name) in mat2_rownames.iter().enumerate() {
        mat2_map.insert(name.as_str(), idx);
    }

    let num_rows = all_rownames.len();
    let num_col1 = mat1.cols();
    let num_col2 = mat2.cols();
    let mut triplets = Vec::with_capacity(mat1.nnz() + mat2.nnz());

    for (out_row, key) in all_rownames.iter().enumerate() {
        if let Some(&src_row) = mat1_map.get(key.as_str()) {
            if let Some(row) = mat1.outer_iterator().nth(src_row) {
                for (col, &val) in row.iter() {
                    triplets.push((out_row, col, val));
                }
            }
        }
        if let Some(&src_row) = mat2_map.get(key.as_str()) {
            if let Some(row) = mat2.outer_iterator().nth(src_row) {
                for (col, &val) in row.iter() {
                    triplets.push((out_row, num_col1 + col, val));
                }
            }
        }
    }

    csc_from_triplets(num_rows, num_col1 + num_col2, &triplets)
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
        let (mean, col_sdev) =
            gene_mean_sdev(&view, &row_index, gene, n_cells, scale, center);
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

pub fn fast_exp_mean_impl(view: CscView<'_>, _display_progress: bool) -> Doubles {
    let nrows = view.nrows as usize;
    let ncols_f = view.ncols as f64;
    let row_index = RowIndex::from_csc_view(&view);

    let rowmeans: Vec<f64> = (0..nrows)
        .map(|row| {
            let sum: f64 = row_index
                .row_range(row)
                .map(|pos| expm1(view.x[row_index.row_x_idx[pos]]))
                .sum();
            log1p(sum / ncols_f)
        })
        .collect();

    Doubles::from_values(rowmeans)
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

pub fn sparse_row_var_std_impl(
    view: CscView<'_>,
    mu: &[f64],
    sd: &[f64],
    vmax: f64,
    _display_progress: bool,
) -> Doubles {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    let row_index = RowIndex::from_csc_view(&view);

    let all_vars: Vec<f64> = (0..n_genes)
        .map(|gene_idx| {
            if sd[gene_idx] == 0.0 {
                return 0.0;
            }
            let range = row_index.row_range(gene_idx);
            let n_zero = n_cells - range.len();
            let mut col_sum = 0.0;
            for pos in range {
                let val = view.x[row_index.row_x_idx[pos]];
                let standardized = ((val - mu[gene_idx]) / sd[gene_idx]).min(vmax);
                col_sum += standardized.powi(2);
            }
            col_sum += ((0.0 - mu[gene_idx]) / sd[gene_idx]).powi(2) * n_zero as f64;
            col_sum / (n_cells as f64 - 1.0)
        })
        .collect();

    Doubles::from_values(all_vars)
}

pub fn fast_log_vmr_impl(view: CscView<'_>, _display_progress: bool) -> Doubles {
    let nrows = view.nrows as usize;
    let ncols = view.ncols as usize;
    let ncols_f = ncols as f64;
    let row_index = RowIndex::from_csc_view(&view);

    let rowdisp: Vec<f64> = (0..nrows)
        .map(|row| {
            let range = row_index.row_range(row);
            let rm: f64 = range
                .clone()
                .map(|pos| expm1(view.x[row_index.row_x_idx[pos]]))
                .sum::<f64>()
                / ncols_f;
            let mut v = 0.0;
            let nn_zero = range.len();
            for pos in range {
                let val = view.x[row_index.row_x_idx[pos]];
                v += (expm1(val) - rm).powi(2);
            }
            v = (v + (ncols - nn_zero) as f64 * rm.powi(2)) / (ncols - 1) as f64;
            (v / rm).ln()
        })
        .collect();

    Doubles::from_values(rowdisp)
}

pub fn sparse_row_var_impl(view: CscView<'_>, _display_progress: bool) -> Doubles {
    let n_genes = view.nrows as usize;
    let n_cells = view.ncols as usize;
    let n_cells_f = n_cells as f64;
    let row_index = RowIndex::from_csc_view(&view);

    let rowdisp: Vec<f64> = (0..n_genes)
        .map(|row| {
            let range = row_index.row_range(row);
            let rm: f64 = range
                .clone()
                .map(|pos| view.x[row_index.row_x_idx[pos]])
                .sum::<f64>()
                / n_cells_f;
            let mut v = 0.0;
            let nn_zero = range.len();
            for pos in range {
                let val = view.x[row_index.row_x_idx[pos]];
                v += (val - rm).powi(2);
            }
            v = (v + (n_cells - nn_zero) as f64 * rm.powi(2)) / (n_cells as f64 - 1.0);
            v
        })
        .collect();

    Doubles::from_values(rowdisp)
}

pub fn replace_cols_impl(
    mat: CscSlots,
    col_idx: &[i32],
    replacement: CscSlots,
) -> CscSlots {
    let nrows = mat.nrows as usize;
    let ncols = mat.ncols as usize;
    let mut triplets = Vec::new();
    let replace_map: HashMap<usize, usize> = col_idx
        .iter()
        .enumerate()
        .map(|(rep_idx, &col)| (col as usize, rep_idx))
        .collect();

    for col in 0..ncols {
        if let Some(&rep_idx) = replace_map.get(&col) {
            for idx in replacement.p[rep_idx] as usize..replacement.p[rep_idx + 1] as usize {
                triplets.push((
                    replacement.i[idx] as usize,
                    col,
                    replacement.x[idx],
                ));
            }
        } else {
            for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
                triplets.push((mat.i[idx] as usize, col, mat.x[idx]));
            }
        }
    }

    csc_from_triplets(nrows, ncols, &triplets)
}

pub fn graph_to_neighbor_helper_impl(mat: CscSlots) -> Robj {
    use crate::sparse::rmatrix_from_ndarray;
    use ndarray::Array2;

    let cs = mat.to_cs_mat();
    let transposed = cs.transpose_view();
    let n_neighbors = transposed
        .outer_iterator()
        .next()
        .map(|row| row.nnz())
        .unwrap_or(0);

    let nrows = transposed.rows();
    let mut nn_idx = Array2::zeros((nrows, n_neighbors));
    let mut nn_dist = Array2::zeros((nrows, n_neighbors));

    for (k, row) in transposed.outer_iterator().enumerate() {
        if row.nnz() != n_neighbors {
            panic!("Not all cells have an equal number of neighbors.");
        }

        let row_idx: Vec<f64> = row.indices().iter().map(|&i| (i + 1) as f64).collect();
        let row_dist: Vec<f64> = row.data().to_vec();

        let mut order: Vec<usize> = (0..row_dist.len()).collect();
        order.sort_by(|&a, &b| {
            row_dist[a]
                .partial_cmp(&row_dist[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        for (i, &ord) in order.iter().enumerate() {
            nn_idx[[k, i]] = row_idx[ord];
            nn_dist[[k, i]] = row_dist[ord];
        }
    }

    Robj::from(vec![
        Robj::from(rmatrix_from_ndarray(nn_idx.view())),
        Robj::from(rmatrix_from_ndarray(nn_dist.view())),
    ])
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
    fn fast_sparse_row_scale_computes_finite_stats() {
        let mat = toy_csc();
        let view = toy_view(&mat);
        let row_index = RowIndex::from_csc_view(&view);
        let (mean, sdev) = gene_mean_sdev(&view, &row_index, 0, 3, true, true);
        assert!(mean.is_finite() && sdev.is_finite() && sdev > 0.0);
    }
}
