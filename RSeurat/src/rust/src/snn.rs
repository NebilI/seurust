use crate::sparse::{csc_slots_from_triplets, dgcmatrix_from_buffers, ndarray_from_rmatrix, CscSlots};
use crate::utils::{row_euclidean_dist, sort_indexes};
use extendr_api::prelude::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicU32, Ordering};

#[cfg(snn_eigen)]
extern "C" {
    fn compute_snn_csc(
        nn_ranked: *const f64,
        nrows: i32,
        ncols: i32,
        prune: f64,
        out_x: *mut *mut f64,
        out_i: *mut *mut i32,
        out_p: *mut *mut i32,
        out_nnz: *mut i32,
        error_msg: *mut std::ffi::c_char,
        error_msg_len: i32,
    ) -> i32;
    fn compute_snn_csc_free(x: *mut f64, i: *mut i32, p: *mut i32);
}

#[cfg(snn_eigen)]
fn compute_snn_eigen_to_r(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    let data = nn_ranked
        .as_robj()
        .as_real_slice()
        .expect("numeric nn_ranked");
    let nrows = nn_ranked.nrows() as i32;
    let ncols = nn_ranked.ncols() as i32;

    let mut out_x: *mut f64 = std::ptr::null_mut();
    let mut out_i: *mut i32 = std::ptr::null_mut();
    let mut out_p: *mut i32 = std::ptr::null_mut();
    let mut out_nnz = 0i32;
    let mut err_buf = vec![0u8; 512];

    let rc = unsafe {
        compute_snn_csc(
            data.as_ptr(),
            nrows,
            ncols,
            prune,
            &mut out_x,
            &mut out_i,
            &mut out_p,
            &mut out_nnz,
            err_buf.as_mut_ptr() as *mut std::ffi::c_char,
            err_buf.len() as i32,
        )
    };

    if rc != 0 {
        let msg = err_buf
            .split(|&b| b == 0)
            .next()
            .unwrap_or(b"compute_snn_csc failed");
        return Err(extendr_api::Error::Other(String::from_utf8_lossy(msg).into_owned()));
    }

    let nnz = out_nnz as usize;
    let n = nrows as usize;
    let mut x_out = Doubles::new(nnz);
    let mut i_out = Integers::new(nnz);
    let mut p_out = Integers::new(n + 1);

    if nnz > 0 {
        x_out
            .as_robj_mut()
            .as_real_slice_mut()
            .expect("numeric x")
            .copy_from_slice(unsafe { std::slice::from_raw_parts(out_x, nnz) });
        i_out
            .as_robj_mut()
            .as_integer_slice_mut()
            .expect("integer i")
            .copy_from_slice(unsafe { std::slice::from_raw_parts(out_i, nnz) });
    }
    p_out
        .as_robj_mut()
        .as_integer_slice_mut()
        .expect("integer p")
        .copy_from_slice(unsafe { std::slice::from_raw_parts(out_p, n + 1) });
    unsafe { compute_snn_csc_free(out_x, out_i, out_p) };

    let dim = Integers::from_values(vec![nrows, nrows]);
    dgcmatrix_from_buffers(x_out, i_out, p_out, dim)
}

fn scale_and_prune(val: f64, k_f: f64, prune: f64) -> Option<f64> {
    let scaled = val / (k_f + (k_f - val));
    if scaled >= prune {
        Some(scaled)
    } else {
        None
    }
}

/// Build reverse-neighbor lists: for each neighbor index v, cells that list v as a neighbor.
fn build_reverse_neighbors(data: &[f64], n_cells: usize, k: usize) -> Vec<Vec<usize>> {
    let mut reverse_neighbors: Vec<Vec<usize>> = vec![Vec::with_capacity(k); n_cells];
    for j in 0..k {
        let base = j * n_cells;
        for i in 0..n_cells {
            let neighbor = data[base + i] as usize - 1;
            reverse_neighbors[neighbor].push(i);
        }
    }
    reverse_neighbors
}

/// Dense shared-neighbor counts with u32 accumulation.
fn count_snn_dense(reverse_neighbors: &[Vec<usize>], n_cells: usize) -> Vec<u32> {
    if n_cells >= 1024 {
        let counts: Vec<AtomicU32> = (0..n_cells * n_cells)
            .map(|_| AtomicU32::new(0))
            .collect();

        reverse_neighbors.par_iter().for_each(|cells| {
            for &i in cells {
                let base = i * n_cells;
                for &j in cells {
                    counts[base + j].fetch_add(1, Ordering::Relaxed);
                }
            }
        });

        return counts
            .into_iter()
            .map(|v| v.into_inner())
            .collect();
    }

    let mut counts = vec![0u32; n_cells * n_cells];
    for cells in reverse_neighbors {
        for &i in cells {
            let base = i * n_cells;
            for &j in cells {
                counts[base + j] += 1;
            }
        }
    }
    counts
}

fn dense_counts_to_triplets(
    counts: &[u32],
    n_cells: usize,
    k_f: f64,
    prune: f64,
) -> Vec<(usize, usize, f64)> {
    let mut triplets = Vec::new();
    for col in 0..n_cells {
        for row in 0..n_cells {
            let val = counts[row * n_cells + col];
            if val > 0 {
                if let Some(scaled) = scale_and_prune(val as f64, k_f, prune) {
                    triplets.push((row, col, scaled));
                }
            }
        }
    }
    triplets
}

fn dense_counts_to_r_dgcmatrix(
    counts: &[u32],
    n_cells: i32,
    k_f: f64,
    prune: f64,
) -> extendr_api::Result<Robj> {
    let n = n_cells as usize;
    let dim = Integers::from_values(vec![n_cells, n_cells]);

    let mut x_vec = Vec::new();
    let mut i_vec = Vec::new();
    let mut p_vec = vec![0i32; n + 1];

    for col in 0..n {
        p_vec[col] = i_vec.len() as i32;
        for row in 0..n {
            let val = counts[row * n + col];
            if val > 0 {
                if let Some(scaled) = scale_and_prune(val as f64, k_f, prune) {
                    i_vec.push(row as i32);
                    x_vec.push(scaled);
                }
            }
        }
    }
    p_vec[n] = i_vec.len() as i32;

    let mut x_out = Doubles::new(x_vec.len());
    let mut i_out = Integers::new(i_vec.len());
    let mut p_out = Integers::new(p_vec.len());

    if !x_vec.is_empty() {
        x_out
            .as_robj_mut()
            .as_real_slice_mut()
            .expect("numeric x")
            .copy_from_slice(&x_vec);
        i_out
            .as_robj_mut()
            .as_integer_slice_mut()
            .expect("integer i")
            .copy_from_slice(&i_vec);
    }
    p_out
        .as_robj_mut()
        .as_integer_slice_mut()
        .expect("integer p")
        .copy_from_slice(&p_vec);

    dgcmatrix_from_buffers(x_out, i_out, p_out, dim)
}

fn triplets_to_csc(n_cells: i32, triplets: Vec<(usize, usize, f64)>) -> CscSlots {
    let n = n_cells as usize;
    let nnz = triplets.len();
    let mut x = Vec::with_capacity(nnz);
    let mut i = Vec::with_capacity(nnz);
    let mut p = vec![0i32; n + 1];
    let mut nz = 0usize;

    for col in 0..n {
        p[col] = nz as i32;
        while nz < nnz && triplets[nz].1 == col {
            let (row, _, val) = triplets[nz];
            i.push(row as i32);
            x.push(val);
            nz += 1;
        }
    }
    p[n] = nz as i32;

    CscSlots {
        x,
        i,
        p,
        nrows: n_cells,
        ncols: n_cells,
    }
}

fn nn_ranked_data(nn_ranked: &RMatrix<f64>) -> &[f64] {
    nn_ranked
        .as_robj()
        .as_real_slice()
        .expect("numeric nn_ranked")
}

/// Core counting kernel: reverse-neighbor lists → sort-reduce → scaled CSC triplets.
pub fn compute_snn_counting_triplets(
    nn_ranked: &RMatrix<f64>,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked_data(nn_ranked);
    let k_f = k as f64;

    let reverse = build_reverse_neighbors(data, n_cells, k);
    let dense = count_snn_dense(&reverse, n_cells);
    let triplets = dense_counts_to_triplets(&dense, n_cells, k_f, prune);
    (n_cells as i32, triplets)
}

/// Compute SNN and return a dgCMatrix with slots written directly in R memory.
pub fn compute_snn_to_r_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    #[cfg(snn_eigen)]
    {
        return compute_snn_eigen_to_r(nn_ranked, prune);
    }

    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked_data(nn_ranked);
    let k_f = k as f64;

    let reverse = build_reverse_neighbors(data, n_cells, k);
    let dense = count_snn_dense(&reverse, n_cells);
    dense_counts_to_r_dgcmatrix(&dense, n_cells as i32, k_f, prune)
}

/// Compute SNN = (neighbor_matrix * neighbor_matrix^T), scaled and pruned.
pub fn compute_snn_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let (n_cells, triplets) = compute_snn_counting_triplets(nn_ranked, prune);
    triplets_to_csc(n_cells, triplets)
}

pub fn write_edge_file_impl(snn: &CscSlots, filename: &str, _display_progress: bool) {
    let mut file = File::create(filename).expect("failed to create edge file");
    let ncols = snn.ncols as usize;
    for col in 0..ncols {
        for idx in snn.p[col] as usize..snn.p[col + 1] as usize {
            let row = snn.i[idx] as usize;
            let val = snn.x[idx];
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
    let mat_arr = ndarray_from_rmatrix(mat);
    let ncols = snn.ncols as usize;
    let mut results = Vec::with_capacity(ncols);

    for col in 0..ncols {
        let start = snn.p[col] as usize;
        let end = snn.p[col + 1] as usize;
        let mut nonzero = Vec::with_capacity(end - start);
        let mut nonzero_idx = Vec::with_capacity(end - start);
        for idx in start..end {
            nonzero.push(snn.x[idx]);
            nonzero_idx.push(snn.i[idx] as usize);
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
                let mut res = row_euclidean_dist(&mat_arr, cell, col);
                if nearest_dist[col] > 0.0 {
                    res -= nearest_dist[col];
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
    use crate::utils::csc_to_dense;

    #[test]
    fn scale_and_prune_matches_formula() {
        let k_f = 20.0;
        let val = 5.0;
        let scaled = scale_and_prune(val, k_f, 0.0).unwrap();
        assert!((scaled - val / (k_f + (k_f - val))).abs() < 1e-12);
    }

    #[test]
    fn counting_kernel_matches_sprs_reference() {
        // Small fixture: 4 cells, k=2 neighbors (1-based indices in columns).
        let n = 4usize;
        let k = 2usize;
        let mut data = vec![0.0; n * k];
        // col 0: cell0->1, cell1->2, cell2->1, cell3->3
        data[0] = 1.0;
        data[1] = 2.0;
        data[2] = 1.0;
        data[3] = 3.0;
        // col 1: cell0->2, cell1->1, cell2->3, cell3->4
        data[4] = 2.0;
        data[5] = 1.0;
        data[6] = 3.0;
        data[7] = 4.0;

        let reverse = build_reverse_neighbors(&data, n, k);
        let dense = count_snn_dense(&reverse, n);
        let triplets = dense_counts_to_triplets(&dense, n, k as f64, 0.0);
        let csc = triplets_to_csc(n as i32, triplets);

        // Reference via sprs SpGEMM on the same neighbor matrix.
        let mut tri = sprs::TriMat::new((n, n));
        for j in 0..k {
            let base = j * n;
            for i in 0..n {
                tri.add_triplet(i, data[base + i] as usize - 1, 1.0);
            }
        }
        let neighbor = tri.to_csc();
        let neighbor_t = neighbor.transpose_view().to_csc();
        let snn = &neighbor * &neighbor_t;

        let mut ref_triplets = Vec::new();
        let k_f = k as f64;
        for (col, col_vec) in snn.outer_iterator().enumerate() {
            for (row, &val) in col_vec.iter() {
                if let Some(scaled) = scale_and_prune(val, k_f, 0.0) {
                    ref_triplets.push((row, col, scaled));
                }
            }
        }
        ref_triplets.sort_by_key(|&(r, c, _)| (c, r));
        let ref_csc = csc_slots_from_triplets(n as i32, n as i32, ref_triplets);

        let dense_got = csc_to_dense(&csc);
        let dense_ref = csc_to_dense(&ref_csc);
        assert_eq!(dense_got.dim(), dense_ref.dim());
        for r in 0..n {
            for c in 0..n {
                assert!(
                    (dense_got[[r, c]] - dense_ref[[r, c]]).abs() < 1e-10,
                    "mismatch at ({r},{c}): {} vs {}",
                    dense_got[[r, c]],
                    dense_ref[[r, c]]
                );
            }
        }
    }
}
