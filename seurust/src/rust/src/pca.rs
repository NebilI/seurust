//! Approximate PCA via the Gram matrix X X' (Seurat 5.6 `EigenGramPCA`).

use crate::sparse::rmatrix_with_data;
use extendr_api::prelude::*;
use nalgebra::{DMatrix, SymmetricEigen};
use rayon::prelude::*;

fn gram_xx_t(x: &[f64], nfeatures: usize, ncells: usize) -> Vec<f64> {
    let mut gram = vec![0.0; nfeatures * nfeatures];
    if nfeatures == 0 || ncells == 0 {
        return gram;
    }
    // gram = X * X^T, X is nfeatures x ncells column-major.
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
    // embeddings = X^T * U, X is nfeatures x ncells, U is nfeatures x npcs.
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
            col.par_iter_mut().for_each(|v| *v = -*v);
        }
    }
}

pub fn eigen_gram_pca_impl(
    object: &RMatrix<f64>,
    npcs: i32,
    weight_by_var: bool,
    _nthreads: i32,
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
    let gram = gram_xx_t(x, nfeatures, ncells);
    let gram_mat = DMatrix::from_column_slice(nfeatures, nfeatures, &gram);
    let eigen = SymmetricEigen::new(gram_mat);

    // nalgebra does not guarantee sorted eigenvalues; Spectra returns them descending.
    let mut order: Vec<usize> = (0..nfeatures).collect();
    order.sort_by(|&a, &b| {
        eigen.eigenvalues[b]
            .partial_cmp(&eigen.eigenvalues[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut d = vec![0.0; npcs];
    let loadings = rmatrix_with_data(nfeatures, npcs, |loadings| {
        for j in 0..npcs {
            let src = order[j];
            d[j] = eigen.eigenvalues[src].max(0.0).sqrt();
            for i in 0..nfeatures {
                loadings[i + j * nfeatures] = eigen.eigenvectors[(i, src)];
            }
        }
        align_loading_signs(loadings, nfeatures, npcs);
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
        // 2 features x 2 cells: columns [1,0] and [0,1] -> I_2, Gram = I_2
        let x = vec![1.0, 0.0, 0.0, 1.0];
        let g = gram_xx_t(&x, 2, 2);
        assert!((g[0] - 1.0).abs() < 1e-12);
        assert!(g[1].abs() < 1e-12);
        assert!(g[2].abs() < 1e-12);
        assert!((g[3] - 1.0).abs() < 1e-12);
    }
}
