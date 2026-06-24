use crate::sparse::{csc_slots_from_sorted_triplets, dgcmatrix_from_buffers, CscSlots};
use extendr_api::prelude::*;
use extendr_ffi::SEXP;
use rayon::prelude::*;
use sprs::{CsMat, TriMat};

#[cfg(snn_eigen)]
#[allow(improper_ctypes)]
extern "C" {
    fn compute_snn_rcpp(nn_ranked: SEXP, prune: f64) -> SEXP;
    fn compute_snn_rcpp_fast(nn_ranked: *const f64, nrows: i32, ncols: i32, prune: f64) -> SEXP;
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
        x_out.as_robj_mut().as_real_slice_mut().expect("numeric x")
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

fn validate_nn_ranked_data(data: &[f64], n_cells: usize, k: usize) -> extendr_api::Result<()> {
    if n_cells == 0 || k == 0 {
        return Err(extendr_api::Error::Other(
            "nn_ranked must have at least one row and one column.".into(),
        ));
    }
    for &idx in data {
        if !idx.is_finite() || idx < 1.0 || idx > n_cells as f64 {
            return Err(extendr_api::Error::Other(format!(
                "nn_ranked contains an invalid neighbor index: {idx}."
            )));
        }
    }
    Ok(())
}

/// Per-cell neighbor lists with duplicate ranks summed, matching Eigen setFromTriplets.
fn neighbors_per_cell_with_counts(
    data: &[f64],
    n_cells: usize,
    k: usize,
) -> Vec<Vec<(usize, u32)>> {
    let mut neighbors = vec![Vec::with_capacity(k); n_cells];
    for rank in 0..k {
        let base = rank * n_cells;
        for i in 0..n_cells {
            let m = data[base + i] as usize - 1;
            let row = &mut neighbors[i];
            if let Some((_, count)) = row.iter_mut().find(|(idx, _)| *idx == m) {
                *count += 1;
            } else {
                row.push((m, 1));
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

    let mut triplets = Vec::with_capacity(col_buckets.iter().map(|entries| entries.len()).sum());
    for (col, mut entries) in col_buckets.into_iter().enumerate() {
        entries.sort_unstable_by_key(|&(row, _)| row);
        for (row, val) in entries {
            triplets.push((row, col, val));
        }
    }

    (n_cells as i32, triplets)
}

/// Build the k-NN neighbor matrix in CSC form, summing duplicate ranks like Eigen.
fn neighbor_csc_from_data(data: &[f64], n_cells: usize, k: usize) -> CsMat<f64> {
    let mut tri = TriMat::new((n_cells, n_cells));
    for rank in 0..k {
        let base = rank * n_cells;
        for i in 0..n_cells {
            tri.add_triplet(i, data[base + i] as usize - 1, 1.0);
        }
    }

    tri.to_csc::<usize>()
}

fn compute_snn_sprs_triplets_from_data(
    data: &[f64],
    n_cells: usize,
    k: usize,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    let neighbor = neighbor_csc_from_data(data, n_cells, k);
    let neighbor_t = neighbor.transpose_view().to_csc();
    let snn = (&neighbor * &neighbor_t).to_csc();
    let k_f = k as f64;
    let mut triplets = Vec::with_capacity(snn.nnz());
    for (col, col_vec) in snn.outer_iterator().enumerate() {
        for (row, &val) in col_vec.iter() {
            if let Some(scaled) = scale_and_prune(val, k_f, prune) {
                triplets.push((row, col, scaled));
            }
        }
    }
    (n_cells as i32, triplets)
}

fn csc_slots_to_triplets(slots: &CscSlots) -> Vec<(usize, usize, f64)> {
    let ncols = slots.ncols as usize;
    let mut triplets = Vec::with_capacity(slots.x.len());
    for col in 0..ncols {
        for idx in slots.p[col] as usize..slots.p[col + 1] as usize {
            triplets.push((slots.i[idx] as usize, col, slots.x[idx]));
        }
    }
    triplets
}

/// Custom SNN kernel that avoids materializing and globally sorting every pair key.
fn compute_snn_accumulating_csc_from_data(
    data: &[f64],
    n_cells: usize,
    k: usize,
    prune: f64,
) -> CscSlots {
    let neighbors = neighbors_per_cell_with_counts(data, n_cells, k);
    let mut inv_count = vec![0usize; n_cells];
    for nbrs in &neighbors {
        for &(m, _) in nbrs {
            inv_count[m] += 1;
        }
    }

    let mut inv_ptr = vec![0usize; n_cells + 1];
    for m in 0..n_cells {
        inv_ptr[m + 1] = inv_ptr[m] + inv_count[m];
    }

    let total = inv_ptr[n_cells];
    let mut inv_cells = vec![0usize; total];
    let mut inv_weights = vec![0u32; total];
    let mut inv_next = inv_ptr.clone();
    for (cell, nbrs) in neighbors.iter().enumerate() {
        for &(m, weight) in nbrs {
            let pos = inv_next[m];
            inv_next[m] += 1;
            inv_cells[pos] = cell;
            inv_weights[pos] = weight;
        }
    }

    let k_f = k as f64;
    let mut accum = vec![0u32; n_cells];
    let mut marks = vec![0u32; n_cells];
    let mut stamp = 1u32;
    let mut touched = Vec::new();
    let mut x = Vec::new();
    let mut i = Vec::new();
    let mut p = Vec::with_capacity(n_cells + 1);

    for col in 0..n_cells {
        p.push(x.len() as i32);
        stamp = stamp.wrapping_add(1);
        if stamp == 0 {
            marks.fill(0);
            stamp = 1;
        }

        for &(m, col_weight) in &neighbors[col] {
            for pos in inv_ptr[m]..inv_ptr[m + 1] {
                let row = inv_cells[pos];
                let weight = inv_weights[pos] * col_weight;
                if marks[row] == stamp {
                    accum[row] += weight;
                } else {
                    marks[row] = stamp;
                    accum[row] = weight;
                    touched.push(row);
                }
            }
        }

        touched.sort_unstable();
        for row in touched.drain(..) {
            if let Some(scaled) = scale_and_prune(accum[row] as f64, k_f, prune) {
                i.push(row as i32);
                x.push(scaled);
            }
        }
    }
    p.push(x.len() as i32);

    CscSlots {
        x,
        i,
        p,
        nrows: n_cells as i32,
        ncols: n_cells as i32,
    }
}

fn compute_snn_accumulating_csc_parallel_from_data(
    data: &[f64],
    n_cells: usize,
    k: usize,
    prune: f64,
) -> CscSlots {
    let neighbors = neighbors_per_cell_with_counts(data, n_cells, k);
    let mut inv_count = vec![0usize; n_cells];
    for nbrs in &neighbors {
        for &(m, _) in nbrs {
            inv_count[m] += 1;
        }
    }

    let mut inv_ptr = vec![0usize; n_cells + 1];
    for m in 0..n_cells {
        inv_ptr[m + 1] = inv_ptr[m] + inv_count[m];
    }

    let total = inv_ptr[n_cells];
    let mut inv_cells = vec![0usize; total];
    let mut inv_weights = vec![0u32; total];
    let mut inv_next = inv_ptr.clone();
    for (cell, nbrs) in neighbors.iter().enumerate() {
        for &(m, weight) in nbrs {
            let pos = inv_next[m];
            inv_next[m] += 1;
            inv_cells[pos] = cell;
            inv_weights[pos] = weight;
        }
    }

    struct Scratch {
        accum: Vec<u32>,
        marks: Vec<u32>,
        touched: Vec<usize>,
        stamp: u32,
    }

    let k_f = k as f64;
    let columns: Vec<(Vec<i32>, Vec<f64>)> = (0..n_cells)
        .into_par_iter()
        .map_init(
            || Scratch {
                accum: vec![0u32; n_cells],
                marks: vec![0u32; n_cells],
                touched: Vec::new(),
                stamp: 1,
            },
            |scratch, col| {
                scratch.stamp = scratch.stamp.wrapping_add(1);
                if scratch.stamp == 0 {
                    scratch.marks.fill(0);
                    scratch.stamp = 1;
                }

                for &(m, col_weight) in &neighbors[col] {
                    for pos in inv_ptr[m]..inv_ptr[m + 1] {
                        let row = inv_cells[pos];
                        let weight = inv_weights[pos] * col_weight;
                        if scratch.marks[row] == scratch.stamp {
                            scratch.accum[row] += weight;
                        } else {
                            scratch.marks[row] = scratch.stamp;
                            scratch.accum[row] = weight;
                            scratch.touched.push(row);
                        }
                    }
                }

                scratch.touched.sort_unstable();
                let mut i_col = Vec::with_capacity(scratch.touched.len());
                let mut x_col = Vec::with_capacity(scratch.touched.len());
                for row in scratch.touched.drain(..) {
                    if let Some(scaled) = scale_and_prune(scratch.accum[row] as f64, k_f, prune) {
                        i_col.push(row as i32);
                        x_col.push(scaled);
                    }
                }
                (i_col, x_col)
            },
        )
        .collect();

    let nnz: usize = columns.iter().map(|(_, x_col)| x_col.len()).sum();
    let mut i = Vec::with_capacity(nnz);
    let mut x = Vec::with_capacity(nnz);
    let mut p = Vec::with_capacity(n_cells + 1);
    for (i_col, x_col) in columns {
        p.push(x.len() as i32);
        i.extend(i_col);
        x.extend(x_col);
    }
    p.push(x.len() as i32);

    CscSlots {
        x,
        i,
        p,
        nrows: n_cells as i32,
        ncols: n_cells as i32,
    }
}

fn compute_snn_best_csc_from_data(data: &[f64], n_cells: usize, k: usize, prune: f64) -> CscSlots {
    if n_cells >= 1000 {
        compute_snn_accumulating_csc_parallel_from_data(data, n_cells, k, prune)
    } else {
        compute_snn_accumulating_csc_from_data(data, n_cells, k, prune)
    }
}

fn compute_snn_accumulating_from_data(
    data: &[f64],
    n_cells: usize,
    k: usize,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    let slots = compute_snn_best_csc_from_data(data, n_cells, k, prune);
    (slots.nrows, csc_slots_to_triplets(&slots))
}

fn compute_snn_spgemm_triplets_from_data(
    data: &[f64],
    n_cells: usize,
    k: usize,
    prune: f64,
) -> (i32, Vec<(usize, usize, f64)>) {
    compute_snn_accumulating_from_data(data, n_cells, k, prune)
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
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked_data(nn_ranked);
    validate_nn_ranked_data(data, n_cells, k)?;
    compute_snn_best_csc_from_data(data, n_cells, k, prune).into_r_dgcmatrix()
}

/// Compute SNN = (neighbor_matrix * neighbor_matrix^T), scaled and pruned.
pub fn compute_snn_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked_data(nn_ranked);
    compute_snn_best_csc_from_data(data, n_cells, k, prune)
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
    use crate::sparse::csc_slots_from_triplets;
    use crate::utils::csc_to_dense;

    fn assert_triplets_close(
        got: &[(usize, usize, f64)],
        expected: &[(usize, usize, f64)],
        tolerance: f64,
    ) {
        assert_eq!(got.len(), expected.len(), "triplet lengths differ");
        for (idx, (g, e)) in got.iter().zip(expected.iter()).enumerate() {
            assert_eq!((g.0, g.1), (e.0, e.1), "triplet index mismatch at {idx}");
            if g.2.is_infinite() || e.2.is_infinite() {
                assert_eq!(g.2, e.2, "triplet value mismatch at {idx}");
                continue;
            }
            assert!(
                (g.2 - e.2).abs() <= tolerance,
                "triplet value mismatch at {idx}: {} vs {}",
                g.2,
                e.2
            );
        }
    }

    fn patterned_nn_data(n: usize, k: usize, seed: usize) -> Vec<f64> {
        let mut data = vec![0.0; n * k];
        for rank in 0..k {
            let base = rank * n;
            for cell in 0..n {
                let val = (cell
                    .wrapping_mul(1_103_515_245)
                    .wrapping_add(rank.wrapping_mul(12_345))
                    .wrapping_add(seed))
                    % n;
                data[base + cell] = (val + 1) as f64;
            }
        }
        data
    }

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

    #[test]
    fn accumulating_kernel_matches_sprs_with_duplicate_neighbor_ranks_and_pruning() {
        let n = 5usize;
        let k = 4usize;
        let data = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 1.0, 4.0, 4.0, 2.0, 3.0, 3.0, 5.0, 5.0, 2.0, 2.0,
            3.0, 5.0, 1.0,
        ];

        for prune in [0.0, 0.01, 0.5] {
            let (_, got) = compute_snn_accumulating_from_data(&data, n, k, prune);
            let (_, expected) = compute_snn_sprs_triplets_from_data(&data, n, k, prune);
            assert_triplets_close(&got, &expected, 1e-12);
        }
    }

    #[test]
    fn accumulating_kernel_matches_sprs_for_larger_patterned_inputs() {
        for &(n, k, seed, prune) in &[(50usize, 20usize, 7usize, 0.0), (200, 20, 19, 0.01)] {
            let data = patterned_nn_data(n, k, seed);
            let (_, got) = compute_snn_accumulating_from_data(&data, n, k, prune);
            let (_, expected) = compute_snn_sprs_triplets_from_data(&data, n, k, prune);
            assert_triplets_close(&got, &expected, 1e-12);
        }
    }

    #[test]
    #[ignore]
    fn benchmark_snn_rust_kernels() {
        use std::time::Instant;

        for &(n, k, reps) in &[(500usize, 20usize, 10usize), (2000, 20, 5)] {
            let data = patterned_nn_data(n, k, 11);

            let start = Instant::now();
            let mut sprs_nnz = 0usize;
            for _ in 0..reps {
                let (_, triplets) = compute_snn_sprs_triplets_from_data(&data, n, k, 0.01);
                sprs_nnz += triplets.len();
            }
            let sprs_elapsed = start.elapsed();

            let start = Instant::now();
            let mut counting_nnz = 0usize;
            for _ in 0..reps {
                let (_, triplets) = compute_snn_counting_from_data(&data, n, k, 0.01);
                counting_nnz += triplets.len();
            }
            let counting_elapsed = start.elapsed();

            let start = Instant::now();
            let mut accumulating_nnz = 0usize;
            for _ in 0..reps {
                let (_, triplets) = compute_snn_accumulating_from_data(&data, n, k, 0.01);
                accumulating_nnz += triplets.len();
            }
            let accumulating_elapsed = start.elapsed();

            assert_eq!(sprs_nnz, counting_nnz);
            assert_eq!(sprs_nnz, accumulating_nnz);
            eprintln!(
                "ComputeSNN Rust kernels n={n}, k={k}, reps={reps}: sprs={:?}, counting={:?}, accumulating={:?}",
                sprs_elapsed, counting_elapsed, accumulating_elapsed
            );
        }
    }
}
