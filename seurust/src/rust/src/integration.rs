use crate::sparse::{csc_from_triplets, CscSlots, CscView};
use crate::utils::sort_indexes;
use extendr_api::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

fn matrix_value(data: &[f64], nrows: usize, row: usize, col: usize) -> f64 {
    data[row + col * nrows]
}

fn col_distance_column_major(data: &[f64], nrows: usize, col_a: usize, col_b: usize) -> f64 {
    let mut sum = 0.0;
    for r in 0..nrows {
        let d = data[r + col_a * nrows] - data[r + col_b * nrows];
        sum += d * d;
    }
    sum.sqrt()
}

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
    let dist_data = distances
        .as_robj()
        .as_real_slice()
        .expect("numeric distances");
    let index_data = cell_index
        .as_robj()
        .as_real_slice()
        .expect("numeric cell_index");
    let dist_nrows = distances.nrows();
    let index_nrows = cell_index.nrows();
    let index_ncols = cell_index.ncols();
    let n_rows = integration_matrix_rownames.len();
    let n_cols = cells2.len();
    let sd_term = (2.0 / sd).powi(2);

    let mut rows_by_name: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, rowname) in integration_matrix_rownames.iter().enumerate() {
        rows_by_name.entry(rowname.as_str()).or_default().push(idx);
    }

    let mut cell_map: Vec<Vec<usize>> = Vec::with_capacity(anchor_cells2.len());
    for name in anchor_cells2 {
        cell_map.push(rows_by_name.get(name.as_str()).cloned().unwrap_or_default());
    }

    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for &cell in cells2 {
        let cell = cell as usize;
        let n_idx = index_ncols;
        let mut k = 0usize;
        for i in 0..n_idx {
            if k >= n_idx {
                break;
            }
            let anchor_idx = matrix_value(index_data, index_nrows, cell, i) as usize - 1;
            if let Some(mnn_idx) = cell_map.get(anchor_idx) {
                for &row in mnn_idx {
                    if k >= n_idx {
                        break;
                    }
                    let dist = matrix_value(dist_data, dist_nrows, cell, i);
                    let to_add = 1.0 - (-dist * anchor_score[row] / sd_term).exp();
                    triplets.push((row, cell, to_add));
                    k += 1;
                }
            }
        }
    }

    if min_dist == 0.0 {
        // Eigen setFromTriplets uses last-wins for duplicate (row, col) entries.
        let mut last_wins: HashMap<(usize, usize), f64> = HashMap::new();
        for &(row, col, val) in &triplets {
            last_wins.insert((row, col), val);
        }
        let sparse_triplets: Vec<(usize, usize, f64)> =
            last_wins.into_iter().map(|((r, c), v)| (r, c, v)).collect();
        let mut mat = csc_from_triplets(n_rows, n_cols, &sparse_triplets);
        let col_sums = mat.col_sums();
        for col in 0..n_cols {
            for idx in mat.p[col] as usize..mat.p[col + 1] as usize {
                mat.x[idx] /= col_sums[col];
            }
        }
        mat
    } else {
        let mut dense = vec![0.0; n_rows * n_cols];
        for col in 0..n_cols {
            for row in 0..n_rows {
                dense[row + col * n_rows] = 1.0 - (-min_dist * anchor_score[row] / sd_term).exp();
            }
        }
        for &(row, col, val) in &triplets {
            dense[row + col * n_rows] = val;
        }
        let mut col_sums = vec![0.0; n_cols];
        for col in 0..n_cols {
            for row in 0..n_rows {
                col_sums[col] += dense[row + col * n_rows];
            }
        }

        let mut x = Vec::with_capacity(n_rows * n_cols);
        let mut i = Vec::with_capacity(n_rows * n_cols);
        let mut p = Vec::with_capacity(n_cols + 1);
        for col in 0..n_cols {
            p.push(x.len() as i32);
            for row in 0..n_rows {
                let value = dense[row + col * n_rows] / col_sums[col];
                if value != 0.0 {
                    i.push(row as i32);
                    x.push(value);
                }
            }
        }
        p.push(x.len() as i32);

        CscSlots {
            x,
            i,
            p,
            nrows: n_rows as i32,
            ncols: n_cols as i32,
        }
    }
}

pub fn integrate_data_impl(
    integration_matrix: CscView<'_>,
    weights: CscView<'_>,
    expression_cells2: CscView<'_>,
) -> CscSlots {
    // expression - t(weights) %*% integration_matrix
    // weights is anchors x genes (CSC), integration_matrix is anchors x cells,
    // expression is genes x cells.
    let n_anchors = weights.nrows as usize;
    let n_genes = weights.ncols as usize;
    let n_cells = expression_cells2.ncols as usize;

    // Flat column-of-transpose: anchor k owns gene/weight pairs in gene order.
    // Adding them in increasing anchor order matches Eigen's sparse product.
    let mut counts = vec![0usize; n_anchors];
    for gene in 0..n_genes {
        let start = weights.p[gene] as usize;
        let end = weights.p[gene + 1] as usize;
        for idx in start..end {
            counts[weights.i[idx] as usize] += 1;
        }
    }
    let mut ptr = vec![0usize; n_anchors + 1];
    for anchor in 0..n_anchors {
        ptr[anchor + 1] = ptr[anchor] + counts[anchor];
    }
    let mut gene_of = vec![0usize; ptr[n_anchors]];
    let mut weight_of = vec![0.0; ptr[n_anchors]];
    let mut cursor = ptr.clone();
    for gene in 0..n_genes {
        let start = weights.p[gene] as usize;
        let end = weights.p[gene + 1] as usize;
        for idx in start..end {
            let anchor = weights.i[idx] as usize;
            let dest = cursor[anchor];
            cursor[anchor] = dest + 1;
            gene_of[dest] = gene;
            weight_of[dest] = weights.x[idx];
        }
    }

    let threads = rayon::current_num_threads().max(1);
    let chunk = n_cells.div_ceil(threads);
    let im_i = integration_matrix.i;
    let im_x = integration_matrix.x;
    let im_p = integration_matrix.p;
    let ex_i = expression_cells2.i;
    let ex_x = expression_cells2.x;
    let ex_p = expression_cells2.p;

    let parts: Vec<(Vec<i32>, Vec<f64>, Vec<i32>)> = (0..threads)
        .into_par_iter()
        .map(|thread| {
            let cell0 = thread * chunk;
            let cell1 = (cell0 + chunk).min(n_cells);
            let mut acc = vec![0.0f64; n_genes];
            let mut col_i = Vec::new();
            let mut col_x = Vec::new();
            let mut col_n = Vec::with_capacity(cell1.saturating_sub(cell0));
            let acc_ptr = acc.as_mut_ptr();
            let gene_ptr = gene_of.as_ptr();
            let weight_ptr = weight_of.as_ptr();
            let im_i_ptr = im_i.as_ptr();
            let im_x_ptr = im_x.as_ptr();
            let ex_i_ptr = ex_i.as_ptr();
            let ex_x_ptr = ex_x.as_ptr();
            for cell in cell0..cell1 {
                let before = col_x.len();
                let im_start = im_p[cell] as usize;
                let im_end = im_p[cell + 1] as usize;
                unsafe {
                    for idx in im_start..im_end {
                        let anchor = *im_i_ptr.add(idx) as usize;
                        let im_val = *im_x_ptr.add(idx);
                        let mut pos = *ptr.get_unchecked(anchor);
                        let w_end = *ptr.get_unchecked(anchor + 1);
                        while pos < w_end {
                            let gene = *gene_ptr.add(pos);
                            *acc_ptr.add(gene) += *weight_ptr.add(pos) * im_val;
                            pos += 1;
                        }
                    }
                    let mut e = *ex_p.get_unchecked(cell) as usize;
                    let e_end = *ex_p.get_unchecked(cell + 1) as usize;
                    for gene in 0..n_genes {
                        let mut expr_val = 0.0;
                        if e < e_end && *ex_i_ptr.add(e) as usize == gene {
                            expr_val = *ex_x_ptr.add(e);
                            e += 1;
                        }
                        let slot = acc_ptr.add(gene);
                        let value = expr_val - *slot;
                        *slot = 0.0;
                        if value != 0.0 {
                            col_i.push(gene as i32);
                            col_x.push(value);
                        }
                    }
                }
                col_n.push((col_x.len() - before) as i32);
            }
            (col_i, col_x, col_n)
        })
        .collect();

    let nnz: usize = parts.iter().map(|(_, x, _)| x.len()).sum();
    let mut x = Vec::with_capacity(nnz);
    let mut i = Vec::with_capacity(nnz);
    let mut p = Vec::with_capacity(n_cells + 1);
    p.push(0);
    for (col_i, col_x, col_n) in parts {
        i.extend(col_i);
        x.extend(col_x);
        for n in col_n {
            p.push(p.last().copied().unwrap_or(0) + n);
        }
    }

    CscSlots {
        x,
        i,
        p,
        nrows: n_genes as i32,
        ncols: n_cells as i32,
    }
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
    let pca_data = query_pca
        .as_robj()
        .as_real_slice()
        .expect("numeric query_pca");
    let qd_data = query_dists
        .as_robj()
        .as_real_slice()
        .expect("numeric query_dists");
    let cn_data = corrected_nns
        .as_robj()
        .as_real_slice()
        .expect("numeric corrected_nns");
    let pca_nrows = query_pca.nrows();
    let qd_nrows = query_dists.nrows();
    let qd_ncols = query_dists.ncols();
    let cn_nrows = corrected_nns.nrows();
    let cn_ncols = corrected_nns.ncols();

    let snn_columns: Vec<(Vec<usize>, Vec<f64>)> = cs
        .outer_iterator()
        .map(|col_vec| {
            let mut nonzero_idx = Vec::new();
            let mut nonzero = Vec::new();
            for (row, &val) in col_vec.iter() {
                nonzero_idx.push(row);
                nonzero.push(val);
            }
            (nonzero_idx, nonzero)
        })
        .collect();

    let scores: Vec<f64> = snn_columns
        .par_iter()
        .enumerate()
        .map(|(i, (nonzero_idx, nonzero))| {
            let order = sort_indexes(&nonzero);
            let mut k_snn_i = k_snn as usize;
            if k_snn_i > order.len() {
                k_snn_i = order.len();
            }

            let mut bw_dists = Vec::new();
            for &ord in &order {
                let cell = nonzero_idx[ord];
                if bw_dists.len() < k_snn_i || nonzero[ord] == nonzero[order[k_snn_i - 1]] {
                    bw_dists.push(col_distance_column_major(pca_data, pca_nrows, cell, i));
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
                matrix_value(qd_data, qd_nrows, i, 1)
            } else {
                0.0
            };
            let bw = bw - first_neighbor_dist;

            let mut q_tps = 0.0;
            for j in 0..qd_ncols {
                q_tps +=
                    (-(matrix_value(qd_data, qd_nrows, i, j) - first_neighbor_dist) / bw).exp();
            }
            q_tps /= qd_ncols as f64;

            let mut c_tps = 0.0;
            for j in 0..cn_ncols {
                let nn_cell = matrix_value(cn_data, cn_nrows, i, j) as usize - 1;
                let dist = col_distance_column_major(pca_data, pca_nrows, i, nn_cell)
                    - first_neighbor_dist;
                c_tps += (-dist / bw).exp();
            }
            c_tps /= cn_ncols as f64;

            c_tps / q_tps
        })
        .collect();

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
        let out = integrate_data_impl(im.as_view(), w.as_view(), expr.as_view());
        assert_eq!(out.x.len(), 1);
        assert!((out.x[0] - 1.5).abs() < 1e-10);
    }
}
