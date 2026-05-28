use crate::sparse::{csc_slots_from_triplets, ndarray_from_rmatrix, CscSlots};
use crate::utils::{row_euclidean_dist, sort_indexes};
use extendr_api::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

const SPRS_SNN_MAX_CELLS: usize = 8192;

fn scale_and_prune(val: f64, k_f: f64, prune: f64) -> Option<f64> {
    let scaled = val / (k_f + (k_f - val));
    if scaled >= prune {
        Some(scaled)
    } else {
        None
    }
}

fn scale_and_prune_to_csc(snn: &sprs::CsMat<f64>, k_f: f64, prune: f64) -> CscSlots {
    let (nrows, ncols) = snn.shape();
    let mut x = Vec::new();
    let mut i = Vec::new();
    let mut p = vec![0i32; ncols + 1];
    let mut nz = 0i32;

    for (col, col_vec) in snn.outer_iterator().enumerate() {
        p[col] = nz;
        for (row, &val) in col_vec.iter() {
            if let Some(scaled) = scale_and_prune(val, k_f, prune) {
                i.push(row as i32);
                x.push(scaled);
                nz += 1;
            }
        }
    }
    p[ncols] = nz;

    CscSlots {
        x,
        i,
        p,
        nrows: nrows as i32,
        ncols: ncols as i32,
    }
}

/// Eigen-equivalent: build binary neighbor matrix, compute A * A^T, scale/prune.
fn compute_snn_sprs_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked.data();
    let k_f = k as f64;

    let nnz = n_cells * k;
    let mut triplets = Vec::with_capacity(nnz);
    for j in 0..k {
        let base = j * n_cells;
        for i in 0..n_cells {
            triplets.push((i, data[base + i] as usize - 1, 1.0));
        }
    }
    let neighbor = csc_slots_from_triplets(n_cells as i32, n_cells as i32, triplets).to_cs_mat();
    let neighbor_t = neighbor.transpose_view().to_csc();
    let snn = &neighbor * &neighbor_t;
    scale_and_prune_to_csc(&snn, k_f, prune)
}

fn collect_scaled_triplets_sparse(
    counts: HashMap<(usize, usize), f64>,
    k_f: f64,
    prune: f64,
) -> Vec<(usize, usize, f64)> {
    let mut triplets = Vec::with_capacity(counts.len());
    for ((row, col), val) in counts {
        if let Some(scaled) = scale_and_prune(val, k_f, prune) {
            triplets.push((row, col, scaled));
        }
    }
    triplets
}

/// Column-wise pair counting for very large n where sprs GEMM is costly.
fn compute_snn_counting_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked.data();
    let k_f = k as f64;

    let mut reverse_neighbors: Vec<Vec<usize>> = vec![Vec::with_capacity(k); n_cells];
    for j in 0..k {
        let base = j * n_cells;
        for i in 0..n_cells {
            let neighbor = data[base + i] as usize - 1;
            reverse_neighbors[neighbor].push(i);
        }
    }

    let mut counts: HashMap<(usize, usize), f64> = HashMap::new();
    for cells in &reverse_neighbors {
        for &i in cells {
            for &j in cells {
                *counts.entry((i, j)).or_insert(0.0) += 1.0;
            }
        }
    }
    let triplets = collect_scaled_triplets_sparse(counts, k_f, prune);
    csc_slots_from_triplets(n_cells as i32, n_cells as i32, triplets)
}

/// Compute SNN = (neighbor_matrix * neighbor_matrix^T), scaled and pruned.
pub fn compute_snn_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let n_cells = nn_ranked.nrows();
    if n_cells <= SPRS_SNN_MAX_CELLS {
        compute_snn_sprs_impl(nn_ranked, prune)
    } else {
        compute_snn_counting_impl(nn_ranked, prune)
    }
}

pub fn write_edge_file_impl(snn: &CscSlots, filename: &str, _display_progress: bool) {
    let cs = snn.to_cs_mat();
    let mut file = File::create(filename).expect("failed to create edge file");
    for (col, col_vec) in cs.outer_iterator().enumerate() {
        for (row, &val) in col_vec.iter() {
            if col >= row {
                continue;
            }
            writeln!(file, "{col}\t{row}\t{val:.15}").unwrap();
        }
    }
}

pub fn direct_snn_to_file_impl(
    nn_ranked: &RMatrix<f64>,
    prune: f64,
    display_progress: bool,
    filename: &str,
) -> CscSlots {
    let snn = compute_snn_impl(nn_ranked, prune);
    write_edge_file_impl(&snn, filename, display_progress);
    snn
}

pub fn snn_smallest_nonzero_dist_impl(
    snn: CscSlots,
    mat: &RMatrix<f64>,
    n: i32,
    nearest_dist: &[f64],
) -> Doubles {
    let cs = snn.to_cs_mat();
    let mat_arr = ndarray_from_rmatrix(mat);
    let mut results = Vec::with_capacity(cs.cols());

    for (i, col_vec) in cs.outer_iterator().enumerate() {
        let mut nonzero = Vec::new();
        let mut nonzero_idx = Vec::new();
        for (row, &val) in col_vec.iter() {
            nonzero.push(val);
            nonzero_idx.push(row);
        }

        let order = sort_indexes(&nonzero);
        let mut n_i = n as usize;
        if n_i > order.len() {
            n_i = order.len();
        }

        let mut dists = Vec::new();
        for &ord in &order {
            let cell = nonzero_idx[ord];
            if dists.len() < n_i || nonzero[ord] == nonzero[order[n_i - 1]] {
                let mut res = row_euclidean_dist(&mat_arr, cell, i);
                if nearest_dist[i] > 0.0 {
                    res -= nearest_dist[i];
                    if res < 0.0 {
                        res = 0.0;
                    }
                }
                dists.push(res);
            } else {
                break;
            }
        }

        let avg_dist = if dists.len() > n_i {
            dists.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            dists[..n_i].iter().sum::<f64>() / n_i as f64
        } else if dists.is_empty() {
            0.0
        } else {
            dists.iter().sum::<f64>() / dists.len() as f64
        };

        results.push(avg_dist);
    }

    Doubles::from_values(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_and_prune_matches_formula() {
        let k_f = 20.0;
        let val = 5.0;
        let scaled = scale_and_prune(val, k_f, 0.0).unwrap();
        assert!((scaled - val / (k_f + (k_f - val))).abs() < 1e-12);
    }
}
