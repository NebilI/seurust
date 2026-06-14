use crate::sparse::{
    csc_slots_from_sorted_triplets, dgcmatrix_from_buffers, dgcmatrix_from_merged_triplets,
    dgcmatrix_from_triplets, CscSlots,
};
use extendr_api::prelude::*;
use extendr_ffi::SEXP;
use sprs::CsMat;

#[cfg(snn_eigen)]
extern "C" {
    fn compute_snn_rcpp(nn_ranked: SEXP, prune: f64) -> SEXP;
    fn compute_snn_rcpp_fast(
        nn_ranked: *const f64,
        nrows: i32,
        ncols: i32,
        prune: f64,
    ) -> SEXP;
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
    fn compute_snn_csc_into(
        nn_ranked: *const f64,
        nrows: i32,
        ncols: i32,
        prune: f64,
        out_x: *mut f64,
        out_x_len: i32,
        out_i: *mut i32,
        out_i_len: i32,
        out_p: *mut i32,
        out_p_len: i32,
        out_nnz_required: *mut i32,
        error_msg: *mut std::ffi::c_char,
        error_msg_len: i32,
    ) -> i32;
    fn compute_snn_csc_clear_cache();
}

#[cfg(snn_eigen)]
fn snn_eigen_error(err_buf: &[u8]) -> extendr_api::Error {
    let msg = err_buf
        .split(|&b| b == 0)
        .next()
        .unwrap_or(b"compute_snn_csc_into failed");
    extendr_api::Error::Other(String::from_utf8_lossy(msg).into_owned())
}

#[cfg(snn_eigen)]
fn compute_snn_eigen_malloc_to_r(
    data: &[f64],
    nrows: i32,
    ncols: i32,
    prune: f64,
) -> extendr_api::Result<Robj> {
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
        return Err(snn_eigen_error(&err_buf));
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

#[cfg(snn_eigen)]
fn compute_snn_eigen_into_to_r(
    data: &[f64],
    nrows: i32,
    ncols: i32,
    prune: f64,
) -> extendr_api::Result<Robj> {
    let n = nrows as usize;
    let mut err_buf = vec![0u8; 512];
    let mut nnz_required = 0i32;

    let sizing_rc = unsafe {
        compute_snn_csc_into(
            data.as_ptr(),
            nrows,
            ncols,
            prune,
            std::ptr::null_mut(),
            0,
            std::ptr::null_mut(),
            0,
            std::ptr::null_mut(),
            0,
            &mut nnz_required,
            err_buf.as_mut_ptr() as *mut std::ffi::c_char,
            err_buf.len() as i32,
        )
    };

    if sizing_rc != -3 {
        unsafe { compute_snn_csc_clear_cache() };
        return Err(snn_eigen_error(&err_buf));
    }

    let nnz = nnz_required.max(0) as usize;
    let mut x_out = Doubles::new(nnz);
    let mut i_out = Integers::new(nnz);
    let mut p_out = Integers::new(n + 1);

    let x_slice = if nnz > 0 {
        x_out
            .as_robj_mut()
            .as_real_slice_mut()
            .expect("numeric x")
    } else {
        &mut []
    };
    let i_slice = if nnz > 0 {
        i_out
            .as_robj_mut()
            .as_integer_slice_mut()
            .expect("integer i")
    } else {
        &mut []
    };
    let p_slice = p_out
        .as_robj_mut()
        .as_integer_slice_mut()
        .expect("integer p");

    let fill_rc = unsafe {
        compute_snn_csc_into(
            data.as_ptr(),
            nrows,
            ncols,
            prune,
            x_slice.as_mut_ptr(),
            nnz as i32,
            i_slice.as_mut_ptr(),
            nnz as i32,
            p_slice.as_mut_ptr(),
            (n + 1) as i32,
            std::ptr::null_mut(),
            err_buf.as_mut_ptr() as *mut std::ffi::c_char,
            err_buf.len() as i32,
        )
    };

    if fill_rc < 0 {
        unsafe { compute_snn_csc_clear_cache() };
        return Err(snn_eigen_error(&err_buf));
    }

    let dim = Integers::from_values(vec![nrows, nrows]);
    dgcmatrix_from_buffers(x_out, i_out, p_out, dim)
}

#[cfg(snn_eigen)]
fn compute_snn_eigen_fast_to_r(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    let data = nn_ranked_data(nn_ranked);
    let out = unsafe {
        compute_snn_rcpp_fast(
            data.as_ptr(),
            nn_ranked.nrows() as i32,
            nn_ranked.ncols() as i32,
            prune,
        )
    };
    Ok(unsafe { Robj::from_sexp(out) })
}

#[cfg(snn_eigen)]
fn compute_snn_eigen_to_r(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    let nn_sexp = unsafe { nn_ranked.as_robj().get() };
    let out = unsafe { compute_snn_rcpp(nn_sexp, prune) };
    Ok(unsafe { Robj::from_sexp(out) })
}

fn scale_and_prune(val: f64, k_f: f64, prune: f64) -> Option<f64> {
    let scaled = val / (k_f + (k_f - val));
    if scaled >= prune {
        Some(scaled)
    } else {
        None
    }
}

fn nn_ranked_data(nn_ranked: &RMatrix<f64>) -> &[f64] {
    nn_ranked
        .as_robj()
        .as_real_slice()
        .expect("numeric nn_ranked")
}

/// Per-cell neighbor lists with duplicate ranks collapsed (matches Eigen setFromTriplets).
fn neighbors_per_cell_deduped(data: &[f64], n_cells: usize, k: usize) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::with_capacity(k); n_cells];
    for rank in 0..k {
        let base = rank * n_cells;
        for i in 0..n_cells {
            let m = data[base + i] as usize - 1;
            let row = &mut neighbors[i];
            if !row.contains(&m) {
                row.push(m);
            }
        }
    }
    neighbors
}

/// SNN[i,l] = sum_m N[i,m] * N[l,m] where N counts neighbor multiplicity.
/// Uses inverted-index pair counting in O(n * k^2) instead of general SpGEMM.
fn compute_snn_counting_from_data(
    data: &[f64],
    n_cells: usize,
    k: usize,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    let mut inv_count = vec![0usize; n_cells];
    for rank in 0..k {
        let base = rank * n_cells;
        for i in 0..n_cells {
            inv_count[data[base + i] as usize - 1] += 1;
        }
    }

    let mut inv_ptr = vec![0usize; n_cells + 1];
    for m in 0..n_cells {
        inv_ptr[m + 1] = inv_ptr[m] + inv_count[m];
    }

    let mut inv_cells = vec![0usize; n_cells * k];
    let mut inv_next = inv_ptr.clone();
    for rank in 0..k {
        let base = rank * n_cells;
        for i in 0..n_cells {
            let m = data[base + i] as usize - 1;
            let pos = inv_next[m];
            inv_next[m] += 1;
            inv_cells[pos] = i;
        }
    }

    let k_f = k as f64;
    let total_pairs: usize = (0..n_cells)
        .map(|m| {
            let len = inv_ptr[m + 1] - inv_ptr[m];
            len * len
        })
        .sum();
    let mut pairs = Vec::with_capacity(total_pairs);
    for m in 0..n_cells {
        let start = inv_ptr[m];
        let end = inv_ptr[m + 1];
        for a in start..end {
            let i = inv_cells[a];
            for b in start..end {
                let l = inv_cells[b];
                pairs.push(((i as u64) << 32) | (l as u64));
            }
        }
    }
    if pairs.is_empty() {
        return (n_cells as i32, Vec::new());
    }

    pairs.sort_unstable();

    let mut row_major: Vec<(usize, usize, f64)> = Vec::new();
    let mut idx = 0;
    while idx < pairs.len() {
        let key = pairs[idx];
        let row = (key >> 32) as usize;
        let col = (key & 0xFFFF_FFFF) as usize;
        let mut count = 1u32;
        idx += 1;
        while idx < pairs.len() && pairs[idx] == key {
            count += 1;
            idx += 1;
        }
        if let Some(scaled) = scale_and_prune(count as f64, k_f, prune) {
            row_major.push((row, col, scaled));
        }
    }

    let mut col_buckets = vec![Vec::new(); n_cells];
    for (row, col, val) in row_major {
        col_buckets[col].push((row, val));
    }

    let mut triplets = Vec::with_capacity(
        col_buckets
            .iter()
            .map(|entries| entries.len())
            .sum(),
    );
    for (col, mut entries) in col_buckets.into_iter().enumerate() {
        entries.sort_unstable_by_key(|&(row, _)| row);
        for (row, val) in entries {
            triplets.push((row, col, val));
        }
    }

    (n_cells as i32, triplets)
}

/// Build the k-NN neighbor matrix directly in CSC form (n_cells x n_cells).
fn neighbor_csc_from_data(data: &[f64], n_cells: usize, k: usize) -> CsMat<f64> {
    let neighbors = neighbors_per_cell_deduped(data, n_cells, k);
    let nnz: usize = neighbors.iter().map(|nbrs| nbrs.len()).sum();
    let mut col_counts = vec![0usize; n_cells];
    for nbrs in &neighbors {
        for &m in nbrs {
            col_counts[m] += 1;
        }
    }

    let mut col_ptr = vec![0usize; n_cells + 1];
    for col in 0..n_cells {
        col_ptr[col + 1] = col_ptr[col] + col_counts[col];
    }

    let mut col_next = col_ptr.clone();
    let mut indices = vec![0usize; nnz];
    let data_vals = vec![1.0f64; nnz];
    for (i, nbrs) in neighbors.iter().enumerate() {
        for &m in nbrs {
            let pos = col_next[m];
            col_next[m] += 1;
            indices[pos] = i;
        }
    }

    CsMat::new_csc((n_cells, n_cells), col_ptr, indices, data_vals)
}

fn compute_snn_spgemm_triplets_from_data(
    data: &[f64],
    n_cells: usize,
    k: usize,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    compute_snn_counting_from_data(data, n_cells, k, prune)
}

fn compute_snn_spgemm_triplets(
    nn_ranked: &RMatrix<f64>,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked_data(nn_ranked);
    compute_snn_spgemm_triplets_from_data(data, n_cells, k, prune)
}

fn triplets_to_csc(n_cells: i32, triplets: Vec<(usize, usize, f64)>) -> CscSlots {
    csc_slots_from_sorted_triplets(n_cells, n_cells, triplets)
}

/// Core counting kernel: SpGEMM → scaled CSC triplets.
pub fn compute_snn_counting_triplets(
    nn_ranked: &RMatrix<f64>,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    compute_snn_spgemm_triplets(nn_ranked, prune)
}

/// Compute SNN and return a dgCMatrix with slots written directly in R memory.
pub fn compute_snn_to_r_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    #[cfg(snn_eigen)]
    {
        return compute_snn_eigen_fast_to_r(nn_ranked, prune);
    }

    #[cfg(not(snn_eigen))]
    {
        let (n_cells, triplets) = compute_snn_spgemm_triplets(nn_ranked, prune);
        dgcmatrix_from_merged_triplets(n_cells, n_cells, triplets)
    }
}

/// Compute SNN = (neighbor_matrix * neighbor_matrix^T), scaled and pruned.
pub fn compute_snn_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let (n_cells, triplets) = compute_snn_spgemm_triplets(nn_ranked, prune);
    triplets_to_csc(n_cells, triplets)
}

pub fn write_edge_file_impl(snn: &CscSlots, filename: &str, _display_progress: bool) {
    use std::fs::File;
    use std::io::Write;

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
    use crate::sparse::ndarray_from_rmatrix;
    use crate::utils::{row_euclidean_dist, sort_indexes};

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
    use crate::sparse::csc_slots_from_triplets;
    use sprs::{CsMat, TriMat};

    #[test]
    fn scale_and_prune_matches_formula() {
        let k_f = 20.0;
        let val = 5.0;
        let scaled = scale_and_prune(val, k_f, 0.0).unwrap();
        assert!((scaled - val / (k_f + (k_f - val))).abs() < 1e-12);
    }

    #[test]
    fn counting_kernel_matches_sprs_with_duplicate_neighbor_ranks() {
        let n = 3usize;
        let k = 3usize;
        // Cell 0 lists neighbor 1 twice; cell 1 lists neighbor 0 once.
        let mut data = vec![0.0; n * k];
        data[0] = 2.0;
        data[1] = 2.0;
        data[2] = 1.0;
        data[3] = 2.0;
        data[4] = 1.0;
        data[5] = 3.0;
        data[6] = 1.0;
        data[7] = 2.0;
        data[8] = 3.0;

        let (_, triplets) = compute_snn_counting_from_data(&data, n, k, 0.0);
        let csc = triplets_to_csc(n as i32, triplets);

        let mut tri = TriMat::new((n, n));
        for j in 0..k {
            let base = j * n;
            for i in 0..n {
                tri.add_triplet(i, data[base + i] as usize - 1, 1.0);
            }
        }
        let neighbor = tri.to_csc::<usize>();
        let neighbor_t = neighbor.transpose_view().to_csc();
        let snn = &neighbor * &neighbor_t;

        let k_f = k as f64;
        let mut ref_triplets = Vec::new();
        for (col, col_vec) in snn.outer_iterator().enumerate() {
            for (row, &val) in col_vec.iter() {
                if let Some(scaled) = scale_and_prune(val, k_f, 0.0) {
                    ref_triplets.push((row, col, scaled));
                }
            }
        }
        ref_triplets.sort_unstable_by_key(|&(r, c, _)| (c, r));
        let ref_csc = csc_slots_from_triplets(n as i32, n as i32, ref_triplets);

        let dense_got = csc_to_dense(&csc);
        let dense_ref = csc_to_dense(&ref_csc);
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

    #[test]
    fn spgemm_kernel_matches_sprs_reference() {
        let n = 4usize;
        let k = 2usize;
        let mut data = vec![0.0; n * k];
        data[0] = 1.0;
        data[1] = 2.0;
        data[2] = 1.0;
        data[3] = 3.0;
        data[4] = 2.0;
        data[5] = 1.0;
        data[6] = 3.0;
        data[7] = 4.0;

        let (_, triplets) = compute_snn_spgemm_triplets_from_data(&data, n, k, 0.0);
        let csc = triplets_to_csc(n as i32, triplets);

        let mut tri = TriMat::new((n, n));
        for j in 0..k {
            let base = j * n;
            for i in 0..n {
                tri.add_triplet(i, data[base + i] as usize - 1, 1.0);
            }
        }
        let neighbor = tri.to_csc::<usize>();
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
        ref_triplets.sort_unstable_by_key(|&(r, c, _)| (c, r));
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
