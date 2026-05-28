use crate::sparse::{csc_from_triplets, rmatrix_from_ndarray, CscSlots, CsrSlots};
use extendr_api::prelude::*;
use extendr_ffi::Rf_runif;
use ndarray::Array2;
use std::collections::HashMap;

fn log1p(x: f64) -> f64 {
    x.ln_1p()
}

fn expm1(x: f64) -> f64 {
    x.exp_m1()
}

pub fn log_norm_impl(mat: &mut CscSlots, scale_factor: i32, _display_progress: bool) {
    let col_sums = mat.col_sums();
    let scale = scale_factor as f64;
    let ncols = mat.ncols as usize;

    for col in 0..ncols {
        for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
            mat.x[idx] = log1p(mat.x[idx] / col_sums[col] * scale);
        }
    }
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

pub fn fast_sparse_row_scale_impl(
    mat: CscSlots,
    scale: bool,
    center: bool,
    scale_max: f64,
    _display_progress: bool,
) -> RMatrix<f64> {
    let n_genes = mat.nrows as usize;
    let n_cells = mat.ncols as usize;
    let transposed = mat.to_cs_mat().transpose_view().to_csc();
    let mut scaled = Array2::zeros((n_genes, n_cells));

    for (gene_idx, col_vec) in transposed.outer_iterator().enumerate() {
        let col_mean: f64 = col_vec.data().iter().sum::<f64>() / n_cells as f64;
        let mut col_sdev = 1.0;

        if scale {
            let mut nn_zero = 0usize;
            let mut var_sum = 0.0;
            if center {
                for &val in col_vec.data() {
                    nn_zero += 1;
                    var_sum += (val - col_mean).powi(2);
                }
                var_sum += col_mean.powi(2) * (n_cells - nn_zero) as f64;
            } else {
                var_sum = col_vec.data().iter().map(|v| v.powi(2)).sum();
            }
            col_sdev = (var_sum / (n_cells - 1) as f64).sqrt();
        }

        let mean = if center { col_mean } else { 0.0 };

        for cell in 0..n_cells {
            let mut value = (0.0 - mean) / col_sdev;
            if value > scale_max {
                value = scale_max;
            }
            scaled[[gene_idx, cell]] = value;
        }
        for (cell, &val) in col_vec.iter() {
            let mut value = (val - mean) / col_sdev;
            if value > scale_max {
                value = scale_max;
            }
            scaled[[gene_idx, cell]] = value;
        }
    }

    rmatrix_from_ndarray(scaled.view())
}

pub fn fast_sparse_row_scale_with_known_stats_impl(
    mat: CscSlots,
    mu: &[f64],
    sigma: &[f64],
    scale: bool,
    center: bool,
    scale_max: f64,
    _display_progress: bool,
) -> RMatrix<f64> {
    let n_genes = mat.nrows as usize;
    let n_cells = mat.ncols as usize;
    let transposed = mat.to_cs_mat().transpose_view().to_csc();
    let mut scaled = Array2::zeros((n_genes, n_cells));

    for (gene_idx, col_vec) in transposed.outer_iterator().enumerate() {
        let col_mean = if center { mu[gene_idx] } else { 0.0 };
        let col_sdev = if scale { sigma[gene_idx] } else { 1.0 };

        for cell in 0..n_cells {
            let mut value = (0.0 - col_mean) / col_sdev;
            if value > scale_max {
                value = scale_max;
            }
            scaled[[gene_idx, cell]] = value;
        }
        for (cell, &val) in col_vec.iter() {
            let mut value = (val - col_mean) / col_sdev;
            if value > scale_max {
                value = scale_max;
            }
            scaled[[gene_idx, cell]] = value;
        }
    }

    rmatrix_from_ndarray(scaled.view())
}

pub fn fast_exp_mean_impl(mat: CscSlots, _display_progress: bool) -> Doubles {
    let nrows = mat.nrows as usize;
    let ncols_f = mat.ncols as f64;
    let cs = mat.to_cs_mat();
    let transposed = cs.transpose_view();
    let mut rowmeans = Vec::with_capacity(nrows);

    for row in transposed.outer_iterator() {
        let rm: f64 = row.data().iter().map(|v| expm1(*v)).sum::<f64>() / ncols_f;
        rowmeans.push(log1p(rm));
    }

    Doubles::from_values(rowmeans)
}

pub fn sparse_row_var2_impl(mat: CscSlots, mu: &[f64], _display_progress: bool) -> Doubles {
    let nrows = mat.nrows as usize;
    let cs = mat.to_cs_mat();
    let transposed = cs.transpose_view();
    let mut all_vars = Vec::with_capacity(transposed.cols());

    for (row_idx, row) in transposed.outer_iterator().enumerate() {
        let mut col_sum = 0.0;
        let mut n_zero = nrows;
        for &val in row.data() {
            n_zero -= 1;
            col_sum += (val - mu[row_idx]).powi(2);
        }
        col_sum += mu[row_idx].powi(2) * n_zero as f64;
        all_vars.push(col_sum / (nrows - 1) as f64);
    }

    Doubles::from_values(all_vars)
}

pub fn sparse_row_var_std_impl(
    mat: CscSlots,
    mu: &[f64],
    sd: &[f64],
    vmax: f64,
    _display_progress: bool,
) -> Doubles {
    let nrows = mat.nrows as usize;
    let ncols = mat.ncols as usize;
    let cs = mat.to_cs_mat();
    let transposed = cs.transpose_view();
    let mut all_vars = vec![0.0; ncols];

    for (row_idx, row) in transposed.outer_iterator().enumerate() {
        if sd[row_idx] == 0.0 {
            continue;
        }
        let mut col_sum = 0.0;
        let mut n_zero = nrows;
        for &val in row.data() {
            n_zero -= 1;
            let standardized = ((val - mu[row_idx]) / sd[row_idx]).min(vmax);
            col_sum += standardized.powi(2);
        }
        col_sum += ((0.0 - mu[row_idx]) / sd[row_idx]).powi(2) * n_zero as f64;
        all_vars[row_idx] = col_sum / (nrows - 1) as f64;
    }

    Doubles::from_values(all_vars)
}

pub fn fast_log_vmr_impl(mat: CscSlots, _display_progress: bool) -> Doubles {
    let nrows = mat.nrows as usize;
    let ncols = mat.ncols as usize;
    let ncols_f = ncols as f64;
    let cs = mat.to_cs_mat();
    let transposed = cs.transpose_view();
    let mut rowdisp = Vec::with_capacity(nrows);

    for row in transposed.outer_iterator() {
        let rm: f64 = row.data().iter().map(|v| expm1(*v)).sum::<f64>() / ncols_f;
        let mut v = 0.0;
        let mut nn_zero = 0usize;
        for &val in row.data() {
            v += (expm1(val) - rm).powi(2);
            nn_zero += 1;
        }
        v = (v + (ncols - nn_zero) as f64 * rm.powi(2)) / (ncols - 1) as f64;
        rowdisp.push((v / rm).ln());
    }

    Doubles::from_values(rowdisp)
}

pub fn sparse_row_var_impl(mat: CscSlots, _display_progress: bool) -> Doubles {
    let nrows = mat.nrows as usize;
    let ncols = mat.ncols as usize;
    let ncols_f = ncols as f64;
    let cs = mat.to_cs_mat();
    let transposed = cs.transpose_view();
    let mut rowdisp = Vec::with_capacity(nrows);

    for row in transposed.outer_iterator() {
        let rm: f64 = row.data().iter().sum::<f64>() / ncols_f;
        let mut v = 0.0;
        let mut nn_zero = 0usize;
        for &val in row.data() {
            v += (val - rm).powi(2);
            nn_zero += 1;
        }
        v = (v + (ncols - nn_zero) as f64 * rm.powi(2)) / (ncols - 1) as f64;
        rowdisp.push(v);
    }

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

    #[test]
    fn log_norm_scales_columns() {
        let mut mat = toy_csc();
        log_norm_impl(&mut mat, 10_000, false);
        assert!(mat.x.iter().all(|v| v.is_finite() && *v >= 0.0));
    }
}
