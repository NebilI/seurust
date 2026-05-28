use crate::sparse::{csc_from_triplets, ndarray_from_rmatrix, CscSlots};
use crate::utils::{col_euclidean_dist, dense_to_csc, sort_indexes};
use extendr_api::prelude::*;
use std::collections::HashMap;

pub fn find_weights_impl(
    cells2: &[i32],
    distances: &RMatrix<f64>,
    anchor_cells2: &[String],
    integration_matrix_rownames: &[String],
    cell_index: &RMatrix<f64>,
    anchor_score: &[f64],
    min_dist: f64,
    sd: f64,
    _display_progress: bool,
) -> CscSlots {
    let dist_mat = ndarray_from_rmatrix(distances);
    let index_mat = ndarray_from_rmatrix(cell_index);
    let n_rows = integration_matrix_rownames.len();
    let n_cols = cells2.len();
    let sd_term = (2.0 / sd).powi(2);

    let mut cell_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, name) in anchor_cells2.iter().enumerate() {
        let matches: Vec<usize> = integration_matrix_rownames
            .iter()
            .enumerate()
            .filter(|(_, rowname)| rowname.as_str() == name.as_str())
            .map(|(idx, _)| idx)
            .collect();
        cell_map.insert(i, matches);
    }

    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for &cell in cells2 {
        let cell = cell as usize;
        let n_idx = index_mat.ncols();
        let mut k = 0usize;
        for i in 0..n_idx {
            if k >= n_idx {
                break;
            }
            let anchor_idx = index_mat[[cell, i]] as usize - 1;
            if let Some(mnn_idx) = cell_map.get(&anchor_idx) {
                for &row in mnn_idx {
                    if k >= n_idx {
                        break;
                    }
                    let dist = dist_mat[[cell, i]];
                    let to_add = 1.0 - (-dist * anchor_score[row] / sd_term).exp();
                    triplets.push((row, cell, to_add));
                    k += 1;
                }
            }
        }
    }

    // Eigen setFromTriplets uses last-wins for duplicate (row, col) entries.
    let mut last_wins: HashMap<(usize, usize), f64> = HashMap::new();
    for &(row, col, val) in &triplets {
        last_wins.insert((row, col), val);
    }
    let triplets: Vec<(usize, usize, f64)> = last_wins
        .into_iter()
        .map(|((r, c), v)| (r, c, v))
        .collect();

    if min_dist == 0.0 {
        let mut mat = csc_from_triplets(n_rows, n_cols, &triplets);
        let col_sums = mat.col_sums();
        for col in 0..n_cols {
            for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
                mat.x[idx] /= col_sums[col];
            }
        }
        mat
    } else {
        let mut dense = ndarray::Array2::<f64>::from_elem((n_rows, n_cols), min_dist);
        for j in 0..n_rows {
            for i in 0..n_cols {
                dense[[j, i]] = 1.0 - (-dense[[j, i]] * anchor_score[j] / sd_term).exp();
            }
        }
        for &(row, col, val) in &triplets {
            dense[[row, col]] = val;
        }
        let col_sums: Vec<f64> = (0..n_cols)
            .map(|c| dense.column(c).sum())
            .collect();
        for col in 0..n_cols {
            for row in 0..n_rows {
                dense[[row, col]] /= col_sums[col];
            }
        }
        dense_to_csc(&dense)
    }
}

pub fn integrate_data_impl(
    integration_matrix: CscSlots,
    weights: CscSlots,
    expression_cells2: CscSlots,
) -> CscSlots {
    let im = integration_matrix.to_cs_mat();
    let w = weights.to_cs_mat();
    let expr = expression_cells2.to_cs_mat();
    let correction = &w.transpose_view().to_csc() * &im;
    let out = &expr - &correction;
    CscSlots::from_cs_mat(&out)
}

pub fn score_helper_impl(
    snn: CscSlots,
    query_pca: &RMatrix<f64>,
    query_dists: &RMatrix<f64>,
    corrected_nns: &RMatrix<f64>,
    k_snn: i32,
    subtract_first_nn: bool,
    _display_progress: bool,
) -> Doubles {
    let cs = snn.to_cs_mat();
    let pca = ndarray_from_rmatrix(query_pca);
    let qd = ndarray_from_rmatrix(query_dists);
    let cn = ndarray_from_rmatrix(corrected_nns);
    let mut scores = Vec::with_capacity(cs.cols());

    for (i, col_vec) in cs.outer_iterator().enumerate() {
        let mut nonzero = Vec::new();
        let mut nonzero_idx = Vec::new();
        for (row, &val) in col_vec.iter() {
            nonzero.push(val);
            nonzero_idx.push(row);
        }

        let order = sort_indexes(&nonzero);
        let mut k_snn_i = k_snn as usize;
        if k_snn_i > order.len() {
            k_snn_i = order.len();
        }

        let mut bw_dists = Vec::new();
        for &ord in &order {
            let cell = nonzero_idx[ord];
            if bw_dists.len() < k_snn_i || nonzero[ord] == nonzero[order[k_snn_i - 1]] {
                bw_dists.push(col_euclidean_dist(&pca, cell, i));
            } else {
                break;
            }
        }

        let bw = if bw_dists.len() > k_snn_i {
            bw_dists.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            bw_dists[..k_snn_i].iter().sum::<f64>() / k_snn_i as f64
        } else if bw_dists.is_empty() {
            0.0
        } else {
            bw_dists.iter().sum::<f64>() / bw_dists.len() as f64
        };

        let first_neighbor_dist = if subtract_first_nn {
            qd[[i, 1]]
        } else {
            0.0
        };
        let bw = bw - first_neighbor_dist;

        let mut q_tps = 0.0;
        for j in 0..qd.ncols() {
            q_tps += (-(qd[[i, j]] - first_neighbor_dist) / bw).exp();
        }
        q_tps /= qd.ncols() as f64;

        let mut c_tps = 0.0;
        for j in 0..cn.ncols() {
            let nn_cell = cn[[i, j]] as usize - 1;
            let dist = col_euclidean_dist(&pca, i, nn_cell) - first_neighbor_dist;
            c_tps += (-dist / bw).exp();
        }
        c_tps /= cn.ncols() as f64;

        scores.push(c_tps / q_tps);
    }

    Doubles::from_values(scores)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::CscSlots;

    #[test]
    fn integrate_data_stays_sparse() {
        let im = CscSlots {
            x: vec![1.0],
            i: vec![0],
            p: vec![0, 1],
            nrows: 1,
            ncols: 1,
        };
        let w = CscSlots {
            x: vec![0.5],
            i: vec![0],
            p: vec![0, 1],
            nrows: 1,
            ncols: 1,
        };
        let expr = CscSlots {
            x: vec![2.0],
            i: vec![0],
            p: vec![0, 1],
            nrows: 1,
            ncols: 1,
        };
        let out = integrate_data_impl(im, w, expr);
        assert_eq!(out.x.len(), 1);
        assert!((out.x[0] - 1.5).abs() < 1e-10);
    }
}
