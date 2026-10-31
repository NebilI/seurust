use crate::sparse::{CscSlots, CscView};
use crate::utils::sort_indexes;
use extendr_api::prelude::*;
use rayon::prelude::*;

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

/// `out = t(left) %*% (right %*% x)` without materializing the cell-by-cell product.
pub fn cca_crossprod_multiply_impl(
    left: &RMatrix<f64>,
    right: &RMatrix<f64>,
    x: &[f64],
) -> extendr_api::Result<Doubles> {
    let nfeatures = left.nrows();
    let ncells = left.ncols();
    if right.nrows() != nfeatures {
        return Err(extendr_api::Error::Other(
            "CCA matrices must have the same number of rows.".into(),
        ));
    }
    if right.ncols() != x.len() {
        return Err(extendr_api::Error::Other(
            "CCA multiplier length must match the number of columns in the right matrix.".into(),
        ));
    }
    let left_data = left.data();
    let right_data = right.data();
    let right_cells = right.ncols();
    let mut feature_scores = vec![0.0; nfeatures];
    for col in 0..right_cells {
        let xv = x[col];
        let base = col * nfeatures;
        for row in 0..nfeatures {
            feature_scores[row] += right_data[base + row] * xv;
        }
    }
    let out: Vec<f64> = (0..ncells)
        .into_par_iter()
        .map(|col| {
            let base = col * nfeatures;
            let mut s = 0.0;
            for row in 0..nfeatures {
                s += left_data[base + row] * feature_scores[row];
            }
            s
        })
        .collect();
    Ok(Doubles::from_values(out))
}

/// Count unique nearest neighbors shared by each anchor pair.
pub fn count_anchor_shared_neighbors_impl(
    indices_aa: &RMatrix<f64>,
    indices_ab: &RMatrix<f64>,
    indices_ba: &RMatrix<f64>,
    indices_bb: &RMatrix<f64>,
    anchor_cell1: &[i32],
    anchor_cell2: &[i32],
    offset: i32,
    k_score: i32,
) -> Integers {
    let nanchors = anchor_cell1.len();
    let max_neighbor = (offset + indices_bb.nrows() as i32) as usize;
    let k_score = k_score.max(0) as usize;
    let aa = indices_aa.data();
    let ab = indices_ab.data();
    let ba = indices_ba.data();
    let bb = indices_bb.data();
    let aa_rows = indices_aa.nrows();
    let ab_rows = indices_ab.nrows();
    let ba_rows = indices_ba.nrows();
    let bb_rows = indices_bb.nrows();
    let offset = offset as usize;

    let scores: Vec<i32> = (0..nanchors)
        .into_par_iter()
        .map(|anchor| {
            let mut neighbors_a = vec![0u32; max_neighbor + 1];
            let mut seen_b = vec![0u32; max_neighbor + 1];
            let stamp = (anchor as u32).wrapping_add(1);
            let anchor_cell_a = (anchor_cell1[anchor] - 1) as usize;
            let anchor_cell_b = (anchor_cell2[anchor] - 1) as usize;
            for k in 0..k_score {
                let aa_id = aa[anchor_cell_a + k * aa_rows] as usize;
                let ab_id = ab[anchor_cell_a + k * ab_rows] as usize + offset;
                if aa_id < neighbors_a.len() {
                    neighbors_a[aa_id] = stamp;
                }
                if ab_id < neighbors_a.len() {
                    neighbors_a[ab_id] = stamp;
                }
            }
            let mut shared = 0i32;
            for k in 0..k_score {
                let neighbor_ba = ba[anchor_cell_b + k * ba_rows] as usize;
                if neighbor_ba < seen_b.len() && seen_b[neighbor_ba] != stamp {
                    seen_b[neighbor_ba] = stamp;
                    if neighbors_a[neighbor_ba] == stamp {
                        shared += 1;
                    }
                }
                let neighbor_bb = bb[anchor_cell_b + k * bb_rows] as usize + offset;
                if neighbor_bb < seen_b.len() && seen_b[neighbor_bb] != stamp {
                    seen_b[neighbor_bb] = stamp;
                    if neighbors_a[neighbor_bb] == stamp {
                        shared += 1;
                    }
                }
            }
            shared
        })
        .collect();
    Integers::from_values(scores)
}

/// Rows of `integration_matrix_rownames` grouped by interned CHARSXP.
struct NameRows {
    keys: Vec<usize>,
    ptr: Vec<u32>,
    rows: Vec<u32>,
}

impl NameRows {
    fn from_names(names: &Strings) -> Self {
        let mut keyed: Vec<(usize, u32)> = names
            .iter()
            .enumerate()
            .map(|(idx, name)| (crate::data_manipulation::char_key(name), idx as u32))
            .collect();
        keyed.sort_unstable_by_key(|&(key, idx)| (key, idx));

        let mut keys = Vec::new();
        let mut ptr = vec![0u32];
        let mut rows = Vec::with_capacity(keyed.len());
        let mut pos = 0usize;
        while pos < keyed.len() {
            let key = keyed[pos].0;
            keys.push(key);
            while pos < keyed.len() && keyed[pos].0 == key {
                rows.push(keyed[pos].1);
                pos += 1;
            }
            ptr.push(rows.len() as u32);
        }
        Self { keys, ptr, rows }
    }

    fn rows_for(&self, key: usize) -> &[u32] {
        match self.keys.binary_search(&key) {
            Ok(uid) => {
                let lo = self.ptr[uid] as usize;
                let hi = self.ptr[uid + 1] as usize;
                &self.rows[lo..hi]
            }
            Err(_) => &[],
        }
    }
}

fn weight_from_dist(dist: f64, score: f64, sd_term: f64) -> f64 {
    1.0 - (-dist * score / sd_term).exp()
}

/// Collect at most `k_max` (row, weight) pairs for one query cell.
fn cell_weight_entries(
    cell: usize,
    dist_data: &[f64],
    dist_nrows: usize,
    index_data: &[f64],
    index_nrows: usize,
    k_max: usize,
    cell_map: &[Vec<u32>],
    anchor_score: &[f64],
    sd_term: f64,
) -> Vec<(u32, f64)> {
    let mut entries = Vec::with_capacity(k_max);
    let mut k = 0usize;
    for neigh in 0..k_max {
        if k >= k_max {
            break;
        }
        let anchor_idx = matrix_value(index_data, index_nrows, cell, neigh) as usize - 1;
        let Some(mnn_idx) = cell_map.get(anchor_idx) else {
            continue;
        };
        let dist = matrix_value(dist_data, dist_nrows, cell, neigh);
        for &row in mnn_idx {
            if k >= k_max {
                break;
            }
            entries.push((
                row,
                weight_from_dist(dist, anchor_score[row as usize], sd_term),
            ));
            k += 1;
        }
    }
    entries
}

/// Eigen `setFromTriplets(..., [](a, b) { return b; })`: keep the last value
/// for each row, then sort rows. `sort_by_key` is stable.
fn last_wins_sorted(entries: &mut Vec<(u32, f64)>) {
    entries.sort_by_key(|&(row, _)| row);
    let mut write = 0usize;
    for read in 0..entries.len() {
        if read + 1 == entries.len() || entries[read].0 != entries[read + 1].0 {
            entries[write] = entries[read];
            write += 1;
        }
    }
    entries.truncate(write);
}

fn normalize_column(entries: &[(u32, f64)]) -> (Vec<i32>, Vec<f64>) {
    let col_sum: f64 = entries.iter().map(|&(_, val)| val).sum();
    let mut i = Vec::with_capacity(entries.len());
    let mut x = Vec::with_capacity(entries.len());
    for &(row, val) in entries {
        i.push(row as i32);
        x.push(val / col_sum);
    }
    (i, x)
}

fn find_weights_from_parts(
    cells2: &[f64],
    dist_data: &[f64],
    dist_nrows: usize,
    index_data: &[f64],
    index_nrows: usize,
    index_ncols: usize,
    cell_map: &[Vec<u32>],
    anchor_score: &[f64],
    n_rows: usize,
    n_cols: usize,
    min_dist: f64,
    sd: f64,
) -> CscSlots {
    let sd_term = (2.0 / sd).powi(2);
    let mut buckets = vec![Vec::<(u32, f64)>::new(); n_cols];
    for &cell_f in cells2 {
        let cell = cell_f as usize;
        if cell >= n_cols {
            continue;
        }
        let extra = cell_weight_entries(
            cell,
            dist_data,
            dist_nrows,
            index_data,
            index_nrows,
            index_ncols,
            cell_map,
            anchor_score,
            sd_term,
        );
        buckets[cell].extend(extra);
    }

    if min_dist == 0.0 {
        let mut i = Vec::new();
        let mut x = Vec::new();
        let mut p = Vec::with_capacity(n_cols + 1);
        for col in 0..n_cols {
            p.push(x.len() as i32);
            last_wins_sorted(&mut buckets[col]);
            let (i_col, x_col) = normalize_column(&buckets[col]);
            i.extend(i_col);
            x.extend(x_col);
        }
        p.push(x.len() as i32);
        CscSlots {
            x,
            i,
            p,
            nrows: n_rows as i32,
            ncols: n_cols as i32,
        }
    } else {
        let mut dense = vec![0.0; n_rows * n_cols];
        for col in 0..n_cols {
            for row in 0..n_rows {
                dense[row + col * n_rows] = weight_from_dist(min_dist, anchor_score[row], sd_term);
            }
        }
        for (col, entries) in buckets.iter().enumerate() {
            for &(row, val) in entries {
                dense[row as usize + col * n_rows] = val;
            }
        }

        let mut x = Vec::with_capacity(n_rows * n_cols);
        let mut i = Vec::with_capacity(n_rows * n_cols);
        let mut p = Vec::with_capacity(n_cols + 1);
        for col in 0..n_cols {
            p.push(x.len() as i32);
            let mut col_sum = 0.0;
            let base = col * n_rows;
            for row in 0..n_rows {
                col_sum += dense[base + row];
            }
            for row in 0..n_rows {
                let value = dense[base + row] / col_sum;
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

pub fn find_weights_impl(
    cells2: &[f64],
    distances: &RMatrix<f64>,
    anchor_cells2: &Strings,
    integration_matrix_rownames: &Strings,
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
    let groups = NameRows::from_names(integration_matrix_rownames);
    let cell_map: Vec<Vec<u32>> = anchor_cells2
        .iter()
        .map(|name| {
            groups
                .rows_for(crate::data_manipulation::char_key(name))
                .to_vec()
        })
        .collect();
    find_weights_from_parts(
        cells2,
        dist_data,
        distances.nrows(),
        index_data,
        cell_index.nrows(),
        cell_index.ncols(),
        &cell_map,
        anchor_score,
        integration_matrix_rownames.len(),
        cells2.len(),
        min_dist,
        sd,
    )
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

    fn csc_to_triplets(slots: &CscSlots) -> Vec<(usize, usize, f64)> {
        let mut out = Vec::with_capacity(slots.x.len());
        for col in 0..slots.ncols as usize {
            for idx in slots.p[col] as usize..slots.p[col + 1] as usize {
                out.push((slots.i[idx] as usize, col, slots.x[idx]));
            }
        }
        out
    }

    #[test]
    fn last_wins_keeps_latest_duplicate_row() {
        let mut entries = vec![(1, 0.1), (0, 0.2), (1, 0.9), (0, 0.4)];
        last_wins_sorted(&mut entries);
        assert_eq!(entries, vec![(0, 0.4), (1, 0.9)]);
    }

    #[test]
    fn find_weights_last_wins_and_normalizes() {
        // One cell, k=3 neighbors. Anchor 0 maps to rows 0 and 1; only the
        // first neighbor is used for both rows, then k hits the cap... wait:
        // k_max=2 so neighbor 0 contributes rows 0 and 1, then k==2.
        let cells2 = [0.0];
        // distances: 1 cell × 2 neighbors, column-major
        let dist_data = [0.5, 1.5];
        // cell_index: neighbor 0 -> anchor 1 (1-based), neighbor 1 unused
        let index_data = [1.0, 2.0];
        // anchor 0 -> rows 0 then 1 (duplicate name group)
        let cell_map = vec![vec![0, 1], vec![1]];
        let anchor_score = [1.0, 0.5];
        let sd = 1.0;
        let sd_term = (2.0_f64 / sd).powi(2);
        let w0 = weight_from_dist(0.5, 1.0, sd_term);
        let w1 = weight_from_dist(0.5, 0.5, sd_term);
        let slots = find_weights_from_parts(
            &cells2,
            &dist_data,
            1,
            &index_data,
            1,
            2,
            &cell_map,
            &anchor_score,
            2,
            1,
            0.0,
            sd,
        );
        let got = csc_to_triplets(&slots);
        let sum = w0 + w1;
        assert_eq!(got.len(), 2);
        assert_eq!((got[0].0, got[0].1), (0, 0));
        assert_eq!((got[1].0, got[1].1), (1, 0));
        assert!((got[0].2 - w0 / sum).abs() < 1e-12);
        assert!((got[1].2 - w1 / sum).abs() < 1e-12);
    }

    #[test]
    fn find_weights_duplicate_row_keeps_last() {
        // Same row emitted twice for one cell; last value wins.
        let cells2 = [0.0];
        let dist_data = [0.2, 0.8];
        // both neighbors point at anchor 0, which maps to row 0 only
        let index_data = [1.0, 1.0];
        let cell_map = vec![vec![0]];
        let anchor_score = [1.0];
        let sd = 1.0;
        let sd_term = (2.0_f64 / sd).powi(2);
        let last = weight_from_dist(0.8, 1.0, sd_term);
        let slots = find_weights_from_parts(
            &cells2,
            &dist_data,
            1,
            &index_data,
            1,
            2,
            &cell_map,
            &anchor_score,
            1,
            1,
            0.0,
            sd,
        );
        let got = csc_to_triplets(&slots);
        assert_eq!(got.len(), 1);
        assert_eq!((got[0].0, got[0].1), (0, 0));
        assert!((got[0].2 - 1.0).abs() < 1e-12);
        assert!(last > 0.0);
    }

    #[test]
    fn find_weights_min_dist_fills_then_overwrites() {
        let cells2 = [0.0];
        let dist_data = [0.4];
        let index_data = [1.0];
        let cell_map = vec![vec![0]];
        let anchor_score = [1.0, 2.0];
        let sd = 1.0;
        let min_dist = 0.25;
        let sd_term = (2.0_f64 / sd).powi(2);
        let fill1 = weight_from_dist(min_dist, 2.0, sd_term);
        let over0 = weight_from_dist(0.4, 1.0, sd_term);
        let slots = find_weights_from_parts(
            &cells2,
            &dist_data,
            1,
            &index_data,
            1,
            1,
            &cell_map,
            &anchor_score,
            2,
            1,
            min_dist,
            sd,
        );
        let got = csc_to_triplets(&slots);
        let sum = over0 + fill1;
        assert_eq!(got.len(), 2);
        assert!((got[0].2 - over0 / sum).abs() < 1e-12);
        assert!((got[1].2 - fill1 / sum).abs() < 1e-12);
    }

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
