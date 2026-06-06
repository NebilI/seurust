use crate::sparse::{csc_slots_from_triplets, dgcmatrix_from_buffers, dgcmatrix_from_triplets, ndarray_from_rmatrix, CscSlots};
use crate::utils::{row_euclidean_dist, sort_indexes};
use extendr_api::prelude::*;
use sprs::TriMat;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

const SPRS_SNN_MAX_CELLS: usize = 8192;

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
struct EigenCsc {
    x: Vec<f64>,
    i: Vec<i32>,
    p: Vec<i32>,
    nrows: i32,
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

#[cfg(snn_eigen)]
fn compute_snn_eigen_csc(nn_ranked: &RMatrix<f64>, prune: f64) -> Result<EigenCsc, String> {
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
        return Err(String::from_utf8_lossy(msg).into_owned());
    }

    let nnz = out_nnz as usize;
    let n = nrows as usize;
    let x = if nnz > 0 {
        unsafe { std::slice::from_raw_parts(out_x, nnz).to_vec() }
    } else {
        Vec::new()
    };
    let i = if nnz > 0 {
        unsafe { std::slice::from_raw_parts(out_i, nnz).to_vec() }
    } else {
        Vec::new()
    };
    let p = unsafe { std::slice::from_raw_parts(out_p, n + 1).to_vec() };
    unsafe { compute_snn_csc_free(out_x, out_i, out_p) };

    Ok(EigenCsc {
        x,
        i,
        p,
        nrows,
    })
}

fn scale_and_prune(val: f64, k_f: f64, prune: f64) -> Option<f64> {
    let scaled = val / (k_f + (k_f - val));
    if scaled >= prune {
        Some(scaled)
    } else {
        None
    }
}

#[cfg(not(snn_eigen))]
fn count_pruned_nnz(snn: &sprs::CsMat<f64>, k_f: f64, prune: f64) -> usize {
    let mut nnz = 0usize;
    for (_, col_vec) in snn.outer_iterator().enumerate() {
        for (_, &val) in col_vec.iter() {
            if scale_and_prune(val, k_f, prune).is_some() {
                nnz += 1;
            }
        }
    }
    nnz
}

#[cfg(not(snn_eigen))]
/// Scale/prune sprs CSC product and write x/i/p directly into preallocated R vectors.
fn scale_and_prune_to_r_dgcmatrix(
    snn: &sprs::CsMat<f64>,
    k_f: f64,
    prune: f64,
) -> extendr_api::Result<Robj> {
    use crate::sparse::dgcmatrix_from_buffers;

    let (nrows, ncols) = snn.shape();
    let nnz = count_pruned_nnz(snn, k_f, prune);
    let dim = Integers::from_values(vec![nrows as i32, ncols as i32]);

    let mut x_out = Doubles::new(nnz);
    let mut i_out = Integers::new(nnz);
    let mut p_out = Integers::new(ncols + 1);

    let x = x_out
        .as_robj_mut()
        .as_real_slice_mut()
        .expect("numeric x");
    let i = i_out
        .as_robj_mut()
        .as_integer_slice_mut()
        .expect("integer i");
    let p = p_out
        .as_robj_mut()
        .as_integer_slice_mut()
        .expect("integer p");

    let mut nz = 0usize;
    for (col, col_vec) in snn.outer_iterator().enumerate() {
        p[col] = nz as i32;
        for (row, &val) in col_vec.iter() {
            if let Some(scaled) = scale_and_prune(val, k_f, prune) {
                i[nz] = row as i32;
                x[nz] = scaled;
                nz += 1;
            }
        }
    }
    p[ncols] = nz as i32;

    dgcmatrix_from_buffers(x_out, i_out, p_out, dim)
}

#[cfg(not(snn_eigen))]
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

#[cfg(not(snn_eigen))]
/// Eigen-equivalent path returning a dgCMatrix with direct slot writes.
fn compute_snn_sprs_to_r(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked.data();
    let k_f = k as f64;

    let mut tri = TriMat::new((n_cells, n_cells));
    for j in 0..k {
        let base = j * n_cells;
        for i in 0..n_cells {
            tri.add_triplet(i, data[base + i] as usize - 1, 1.0);
        }
    }
    let neighbor = tri.to_csc();
    let neighbor_t = neighbor.transpose_view().to_csc();
    let snn = &neighbor * &neighbor_t;
    scale_and_prune_to_r_dgcmatrix(&snn, k_f, prune)
}

#[cfg(not(snn_eigen))]
/// Column-wise pair counting for very large n.
fn compute_snn_counting_to_r(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
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
    dgcmatrix_from_triplets(n_cells as i32, n_cells as i32, triplets)
}

/// Compute SNN and return a dgCMatrix with slots written directly in R memory.
pub fn compute_snn_to_r_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> extendr_api::Result<Robj> {
    #[cfg(snn_eigen)]
    {
        return compute_snn_eigen_to_r(nn_ranked, prune);
    }
    #[cfg(not(snn_eigen))]
    {
        let n_cells = nn_ranked.nrows();
        if n_cells <= SPRS_SNN_MAX_CELLS {
            compute_snn_sprs_to_r(nn_ranked, prune)
        } else {
            compute_snn_counting_to_r(nn_ranked, prune)
        }
    }
}

/// Compute SNN = (neighbor_matrix * neighbor_matrix^T), scaled and pruned.
pub fn compute_snn_impl(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    #[cfg(snn_eigen)]
    {
        let csc = compute_snn_eigen_csc(nn_ranked, prune).expect("compute_snn_eigen_csc");
        return CscSlots {
            x: csc.x,
            i: csc.i,
            p: csc.p,
            nrows: csc.nrows,
            ncols: csc.nrows,
        };
    }

    #[cfg(not(snn_eigen))]
    compute_snn_impl_rust(nn_ranked, prune)
}

#[cfg(not(snn_eigen))]
fn compute_snn_impl_rust(nn_ranked: &RMatrix<f64>, prune: f64) -> CscSlots {
    let n_cells = nn_ranked.nrows();
    let k = nn_ranked.ncols();
    let data = nn_ranked.data();
    let k_f = k as f64;

    if n_cells <= SPRS_SNN_MAX_CELLS {
        let mut tri = TriMat::new((n_cells, n_cells));
        for j in 0..k {
            let base = j * n_cells;
            for i in 0..n_cells {
                tri.add_triplet(i, data[base + i] as usize - 1, 1.0);
            }
        }
        let neighbor = tri.to_csc();
        let neighbor_t = neighbor.transpose_view().to_csc();
        let snn = &neighbor * &neighbor_t;
        scale_and_prune_to_csc(&snn, k_f, prune)
    } else {
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
}

#[cfg(not(snn_eigen))]
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
