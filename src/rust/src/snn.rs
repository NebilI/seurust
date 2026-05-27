use crate::sparse::{csc_from_triplets, ndarray_from_rmatrix, CscSlots};
use crate::utils::{row_euclidean_dist, sort_indexes};
use extendr_api::prelude::*;
use sprs::TriMat;
use std::fs::File;
use std::io::Write;

pub fn compute_snn_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let nn = ndarray_from_rmatrix(nn_ranked);
    let (n_cells, k) = nn.dim();
    let mut triplets = Vec::with_capacity(n_cells * k);

    for j in 0..k {
        for i in 0..n_cells {
            let col = nn[[i, j]] as usize - 1;
            triplets.push((i, col, 1.0));
        }
    }

    let mut tri = TriMat::new((n_cells, n_cells));
    for &(r, c, v) in &triplets {
        tri.add_triplet(r, c, v);
    }
    let snn = tri.to_csc::<usize>();
    let snn_t = snn.transpose_view().to_csc();
    let prod = &snn * &snn_t;

    let k_f = k as f64;
    let mut out_triplets = Vec::new();
    for (col, col_vec) in prod.outer_iterator().enumerate() {
        for (row, &val) in col_vec.iter() {
            let scaled = val / (k_f + (k_f - val));
            if scaled >= prune {
                out_triplets.push((row, col, scaled));
            }
        }
    }

    csc_from_triplets(n_cells, n_cells, &out_triplets)
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
