//! Approximate PCA via the Gram matrix X X' (Seurat 5.6 `EigenGramPCA`).

use crate::sparse::rmatrix_with_data;
use crate::utils::use_parallel;
use extendr_api::prelude::*;
use nalgebra::{DMatrix, SymmetricEigen};
use rayon::prelude::*;

#[derive(Clone, Copy)]
struct SendPtr(*mut f64);
unsafe impl Send for SendPtr {}
unsafe impl Sync for SendPtr {}

fn gram_thread_count(nthreads: i32, nfeatures: usize) -> usize {
    let want = if nthreads <= 0 {
        rayon::current_num_threads()
    } else {
        nthreads as usize
    };
    want.max(1).min(nfeatures.max(1))
}

fn gram_xx_t_serial(x: &[f64], nfeatures: usize, ncells: usize) -> Vec<f64> {
    let mut gram = vec![0.0; nfeatures * nfeatures];
    if nfeatures == 0 || ncells == 0 {
        return gram;
    }
    unsafe {
        matrixmultiply::dgemm(
            nfeatures,
            ncells,
            nfeatures,
            1.0,
            x.as_ptr(),
            1,
            nfeatures as isize,
            x.as_ptr(),
            nfeatures as isize,
            1,
            0.0,
            gram.as_mut_ptr(),
            1,
            nfeatures as isize,
        );
    }
    gram
}

/// Symmetric `X Xᵀ`: lower block triangle in parallel, then mirrored.
/// Each output block keeps the same inner-product order as a full `dgemm`.
fn gram_xx_t_blocked(x: &[f64], nfeatures: usize, ncells: usize, n_blocks: usize) -> Vec<f64> {
    let mut gram = vec![0.0; nfeatures * nfeatures];
    if nfeatures == 0 || ncells == 0 {
        return gram;
    }
    let n_blocks = n_blocks.max(1).min(nfeatures);
    let bs = nfeatures.div_ceil(n_blocks);
    let pairs: Vec<(usize, usize)> = (0..n_blocks)
        .flat_map(|bj| (bj..n_blocks).map(move |bi| (bi, bj)))
        .collect();
    let gp = SendPtr(gram.as_mut_ptr());
    pairs.par_iter().for_each(|&(bi, bj)| {
        let gp = gp;
        let i0 = bi * bs;
        let j0 = bj * bs;
        if i0 >= nfeatures || j0 >= nfeatures {
            return;
        }
        let h = bs.min(nfeatures - i0);
        let w = bs.min(nfeatures - j0);
        unsafe {
            matrixmultiply::dgemm(
                h,
                ncells,
                w,
                1.0,
                x.as_ptr().add(i0),
                1,
                nfeatures as isize,
                x.as_ptr().add(j0),
                nfeatures as isize,
                1,
                0.0,
                gp.0.add(i0 + j0 * nfeatures),
                1,
                nfeatures as isize,
            );
        }
    });
    for j in 0..nfeatures {
        for i in 0..j {
            gram[i + j * nfeatures] = gram[j + i * nfeatures];
        }
    }
    gram
}

fn gram_xx_t(x: &[f64], nfeatures: usize, ncells: usize, nthreads: i32) -> Vec<f64> {
    let work = nfeatures.saturating_mul(ncells).saturating_mul(nfeatures);
    let n_blocks = gram_thread_count(nthreads, nfeatures);
    if n_blocks <= 1 || !use_parallel(n_blocks, work) {
        gram_xx_t_serial(x, nfeatures, ncells)
    } else {
        gram_xx_t_blocked(x, nfeatures, ncells, n_blocks)
    }
}

fn embeddings_xt_u_into(
    x: &[f64],
    loadings: &[f64],
    nfeatures: usize,
    ncells: usize,
    npcs: usize,
    embeddings: &mut [f64],
) {
    if nfeatures == 0 || ncells == 0 || npcs == 0 {
        return;
    }
    debug_assert_eq!(embeddings.len(), ncells * npcs);
    unsafe {
        matrixmultiply::dgemm(
            ncells,
            nfeatures,
            npcs,
            1.0,
            x.as_ptr(),
            nfeatures as isize,
            1,
            loadings.as_ptr(),
            1,
            nfeatures as isize,
            0.0,
            embeddings.as_mut_ptr(),
            1,
            ncells as isize,
        );
    }
}

/// Align each eigenvector so the entry of largest magnitude is non-negative.
/// Spectra and nalgebra can differ by a global sign per column.
fn align_loading_signs(loadings: &mut [f64], nfeatures: usize, npcs: usize) {
    for j in 0..npcs {
        let col = &mut loadings[j * nfeatures..(j + 1) * nfeatures];
        let mut best_idx = 0usize;
        let mut best_abs = 0.0;
        for (i, &v) in col.iter().enumerate() {
            let a = v.abs();
            if a > best_abs {
                best_abs = a;
                best_idx = i;
            }
        }
        if col[best_idx] < 0.0 {
            for v in col.iter_mut() {
                *v = -*v;
            }
        }
    }
}

fn top_eigenpairs_nalgebra(gram: &[f64], nfeatures: usize, npcs: usize) -> (Vec<f64>, Vec<f64>) {
    let gram_mat = DMatrix::from_column_slice(nfeatures, nfeatures, gram);
    let eigen = SymmetricEigen::new(gram_mat);
    let mut order: Vec<usize> = (0..nfeatures).collect();
    order.sort_by(|&a, &b| {
        eigen.eigenvalues[b]
            .partial_cmp(&eigen.eigenvalues[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut values = vec![0.0; npcs];
    let mut vectors = vec![0.0; nfeatures * npcs];
    for (j, &src) in order.iter().take(npcs).enumerate() {
        values[j] = eigen.eigenvalues[src];
        for i in 0..nfeatures {
            vectors[i + j * nfeatures] = eigen.eigenvectors[(i, src)];
        }
    }
    (values, vectors)
}

/// Top `npcs` eigenpairs of a symmetric Gram via LAPACK `dsyevr(RANGE='I')`.
/// Falls back to nalgebra if the solver does not converge.
#[cfg(feature = "r-lapack")]
fn top_eigenpairs(gram: &[f64], nfeatures: usize, npcs: usize) -> (Vec<f64>, Vec<f64>) {
    match top_eigenpairs_dsyevr(gram, nfeatures, npcs) {
        Some(result) => result,
        None => top_eigenpairs_nalgebra(gram, nfeatures, npcs),
    }
}

#[cfg(not(feature = "r-lapack"))]
fn top_eigenpairs(gram: &[f64], nfeatures: usize, npcs: usize) -> (Vec<f64>, Vec<f64>) {
    top_eigenpairs_nalgebra(gram, nfeatures, npcs)
}

#[cfg(feature = "r-lapack")]
fn top_eigenpairs_dsyevr(gram: &[f64], n: usize, npcs: usize) -> Option<(Vec<f64>, Vec<f64>)> {
    use std::os::raw::c_char;

    extern "C" {
        fn dsyevr_(
            jobz: *const c_char,
            range: *const c_char,
            uplo: *const c_char,
            n: *const i32,
            a: *mut f64,
            lda: *const i32,
            vl: *const f64,
            vu: *const f64,
            il: *const i32,
            iu: *const i32,
            abstol: *const f64,
            m: *mut i32,
            w: *mut f64,
            z: *mut f64,
            ldz: *const i32,
            isuppz: *mut i32,
            work: *mut f64,
            lwork: *const i32,
            iwork: *mut i32,
            liwork: *const i32,
            info: *mut i32,
        );
    }

    if n == 0 || npcs == 0 || npcs > n {
        return None;
    }
    let n_i = n as i32;
    let jobz = b"V\0".as_ptr() as *const c_char;
    let range = b"I\0".as_ptr() as *const c_char;
    let uplo = b"U\0".as_ptr() as *const c_char;
    let mut a = gram.to_vec();
    let lda = n_i;
    let vl = 0.0;
    let vu = 0.0;
    let il = (n - npcs + 1) as i32;
    let iu = n_i;
    let abstol = 0.0;
    let mut m = 0i32;
    let mut w = vec![0.0; n];
    let mut z = vec![0.0; n * npcs];
    let ldz = n_i;
    let mut isuppz = vec![0i32; 2 * npcs];
    let mut work0 = [0.0; 1];
    let mut iwork0 = [0i32; 1];
    let mut lwork = -1i32;
    let mut liwork = -1i32;
    let mut info = 0i32;
    unsafe {
        dsyevr_(
            jobz,
            range,
            uplo,
            &n_i,
            a.as_mut_ptr(),
            &lda,
            &vl,
            &vu,
            &il,
            &iu,
            &abstol,
            &mut m,
            w.as_mut_ptr(),
            z.as_mut_ptr(),
            &ldz,
            isuppz.as_mut_ptr(),
            work0.as_mut_ptr(),
            &lwork,
            iwork0.as_mut_ptr(),
            &liwork,
            &mut info,
        );
    }
    if info != 0 {
        return None;
    }
    lwork = work0[0] as i32;
    liwork = iwork0[0];
    if lwork < 1 || liwork < 1 {
        return None;
    }
    let mut work = vec![0.0; lwork as usize];
    let mut iwork = vec![0i32; liwork as usize];
    // `dsyevr` destroys A; query used a copy, so restore.
    a.copy_from_slice(gram);
    m = 0;
    info = 0;
    unsafe {
        dsyevr_(
            jobz,
            range,
            uplo,
            &n_i,
            a.as_mut_ptr(),
            &lda,
            &vl,
            &vu,
            &il,
            &iu,
            &abstol,
            &mut m,
            w.as_mut_ptr(),
            z.as_mut_ptr(),
            &ldz,
            isuppz.as_mut_ptr(),
            work.as_mut_ptr(),
            &lwork,
            iwork.as_mut_ptr(),
            &liwork,
            &mut info,
        );
    }
    if info != 0 || m != npcs as i32 {
        return None;
    }
    // `dsyevr` returns eigenvalues in ascending order; we want descending.
    let mut values = vec![0.0; npcs];
    let mut vectors = vec![0.0; n * npcs];
    for j in 0..npcs {
        let src = npcs - 1 - j;
        values[j] = w[src];
        for i in 0..n {
            vectors[i + j * n] = z[i + src * n];
        }
    }
    Some((values, vectors))
}

pub fn eigen_gram_pca_impl(
    object: &RMatrix<f64>,
    npcs: i32,
    weight_by_var: bool,
    nthreads: i32,
) -> extendr_api::Result<Robj> {
    let nfeatures = object.nrows();
    let ncells = object.ncols();
    if npcs < 1 || (npcs as usize) >= nfeatures.min(ncells) {
        return Err(extendr_api::Error::Other(
            "npcs must be positive and strictly less than both matrix dimensions.".into(),
        ));
    }
    let npcs = npcs as usize;
    let x = object.data();
    let gram = gram_xx_t(x, nfeatures, ncells, nthreads);
    let (values, mut vectors) = top_eigenpairs(&gram, nfeatures, npcs);
    align_loading_signs(&mut vectors, nfeatures, npcs);

    let mut d = vec![0.0; npcs];
    let loadings = rmatrix_with_data(nfeatures, npcs, |loadings| {
        loadings.copy_from_slice(&vectors);
        for j in 0..npcs {
            d[j] = values[j].max(0.0).sqrt();
        }
    });

    let mut embeddings = rmatrix_with_data(ncells, npcs, |embeddings| {
        embeddings_xt_u_into(x, loadings.data(), nfeatures, ncells, npcs, embeddings);
    });
    if !weight_by_var {
        let max_d = d.iter().copied().fold(0.0_f64, f64::max);
        let min_singular = max_d * f64::EPSILON.sqrt();
        if d.iter().copied().fold(f64::INFINITY, f64::min) <= min_singular {
            return Err(extendr_api::Error::Other(
                "Unweighted PCA requires SVD for numerically zero singular values.".into(),
            ));
        }
        let embeddings_data = embeddings.data_mut();
        for j in 0..npcs {
            let scale = d[j];
            for i in 0..ncells {
                embeddings_data[i + j * ncells] /= scale;
            }
        }
    }

    let denom = (ncells.max(2) - 1) as f64;
    let sdev: Vec<f64> = d.iter().map(|v| v / denom.sqrt()).collect();
    Ok(Robj::from(list!(
        loadings = loadings,
        embeddings = embeddings,
        sdev = Doubles::from_values(sdev)
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gram_of_identity_block_is_ones() {
        let x = vec![1.0, 0.0, 0.0, 1.0];
        let g = gram_xx_t(&x, 2, 2, 1);
        assert!((g[0] - 1.0).abs() < 1e-12);
        assert!(g[1].abs() < 1e-12);
        assert!(g[2].abs() < 1e-12);
        assert!((g[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn blocked_gram_matches_serial_bits() {
        let mut x = vec![0.0; 12 * 20];
        for (i, v) in x.iter_mut().enumerate() {
            *v = (i % 7) as f64 * 0.1 - 0.3;
        }
        let serial = gram_xx_t_serial(&x, 12, 20);
        let blocked = gram_xx_t_blocked(&x, 12, 20, 3);
        for (a, b) in serial.iter().zip(&blocked) {
            assert_eq!(a.to_bits(), b.to_bits(), "{a} vs {b}");
        }
    }

    #[test]
    fn nalgebra_top_pairs_match_full_sort() {
        let x: Vec<f64> = (0..8 * 15).map(|i| ((i * 17) % 13) as f64 * 0.05).collect();
        let g = gram_xx_t_serial(&x, 8, 15);
        let (vals, vecs) = top_eigenpairs_nalgebra(&g, 8, 3);
        let full = SymmetricEigen::new(DMatrix::from_column_slice(8, 8, &g));
        let mut order: Vec<usize> = (0..8).collect();
        order.sort_by(|&a, &b| {
            full.eigenvalues[b]
                .partial_cmp(&full.eigenvalues[a])
                .unwrap()
        });
        for j in 0..3 {
            assert!((vals[j] - full.eigenvalues[order[j]]).abs() < 1e-12);
            let mut same = 0.0;
            let mut flipped = 0.0;
            for i in 0..8 {
                let e = full.eigenvectors[(i, order[j])];
                same += (vecs[i + j * 8] - e).abs();
                flipped += (vecs[i + j * 8] + e).abs();
            }
            assert!(same.min(flipped) < 1e-10);
        }
    }
}
