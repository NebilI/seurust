use crate::sparse::rmatrix_from_column_major;
use extendr_api::prelude::*;
use rayon::prelude::*;

pub fn standardize_impl(mat: &RMatrix<f64>, _display_progress: bool) -> RMatrix<f64> {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    let values = mat.data();
    let mut out = vec![0.0; nrows * ncols];

    out.par_chunks_mut(nrows.max(1))
        .enumerate()
        .for_each(|(c, out_col)| {
            let start = c * nrows;
            let col = &values[start..start + nrows];
            let col_mean = if nrows == 0 {
                0.0
            } else {
                col.iter().sum::<f64>() / nrows as f64
            };
            let col_sdev = (col
                .iter()
                .map(|v| {
                    let d = v - col_mean;
                    d * d
                })
                .sum::<f64>()
                / (nrows as f64 - 1.0))
                .sqrt();
            for r in 0..nrows {
                out_col[r] = (col[r] - col_mean) / col_sdev;
            }
        });

    rmatrix_from_column_major(&out, nrows, ncols)
}

/// Column-major copy of `values` (nrows x ncols), optionally with each column
/// centered on its mean (Eigen's `mat.rowwise() - mat.colwise().mean()`).
fn column_major_copy(values: &[f64], nrows: usize, ncols: usize, center: bool) -> Vec<f64> {
    let mut out = values.to_vec();
    if center && nrows > 0 {
        out.par_chunks_mut(nrows).for_each(|col| {
            let mean = col.iter().sum::<f64>() / nrows as f64;
            for v in col.iter_mut() {
                *v -= mean;
            }
        });
    }
    debug_assert_eq!(out.len(), nrows * ncols);
    out
}

/// Width of the output column blocks handed to each rayon task.
fn gemm_block_width(ncols_out: usize) -> usize {
    let tasks = rayon::current_num_threads().max(1) * 4;
    ncols_out.div_ceil(tasks).max(32)
}

/// out (m x ncols_b, column-major, leading dim m) = a^T * b[:, cols], where
/// `a` is n x m and `b` is n x ncols_b, both column-major. Only the first
/// `rows_needed(block_start, block_end)` rows of each output block are computed.
fn crossprod_blocks(
    a: &[f64],
    b: &[f64],
    n: usize,
    m: usize,
    ncols_b: usize,
    out: &mut [f64],
    rows_needed: impl Fn(usize, usize) -> usize + Sync,
) {
    if m == 0 || ncols_b == 0 {
        return;
    }
    let block = gemm_block_width(ncols_b);
    out.par_chunks_mut(m * block)
        .enumerate()
        .for_each(|(blk, out_block)| {
            let j0 = blk * block;
            let width = out_block.len() / m;
            let rows = rows_needed(j0, j0 + width).min(m);
            if n == 0 {
                out_block.fill(0.0);
                return;
            }
            // SAFETY: all pointers and strides describe in-bounds column-major
            // buffers: a^T is (rows x n) with row stride n, b[:, j0..] is
            // (n x width) with column stride n, and out_block is (rows x width)
            // with column stride m.
            unsafe {
                matrixmultiply::dgemm(
                    rows,
                    n,
                    width,
                    1.0,
                    a.as_ptr(),
                    n as isize,
                    1,
                    b.as_ptr().add(j0 * n),
                    1,
                    n as isize,
                    0.0,
                    out_block.as_mut_ptr(),
                    1,
                    m as isize,
                );
            }
        });
}

pub fn fast_cov_impl(mat: &RMatrix<f64>, center: bool) -> RMatrix<f64> {
    let n = mat.nrows();
    let p = mat.ncols();
    let values = column_major_copy(mat.data(), n, p, center);
    let denom = n as f64 - 1.0;

    let mut out = RMatrix::<f64>::new(p, p);
    if p == 0 {
        return out;
    }
    let cov = out.data_mut();
    // X^T X is symmetric: compute the upper triangle block-wise, then mirror.
    crossprod_blocks(&values, &values, n, p, p, cov, |_, j1| j1);
    for j in 0..p {
        for i in 0..=j {
            let v = cov[i + j * p] / denom;
            cov[i + j * p] = v;
            cov[j + i * p] = v;
        }
    }
    out
}

pub fn fast_cov_mats_impl(
    mat1: &RMatrix<f64>,
    mat2: &RMatrix<f64>,
    center: bool,
) -> extendr_api::Result<RMatrix<f64>> {
    let n = mat1.nrows();
    if mat2.nrows() != n {
        return Err(extendr_api::Error::Other(format!(
            "FastCovMats: mat1 has {n} rows but mat2 has {}.",
            mat2.nrows()
        )));
    }
    let p1 = mat1.ncols();
    let p2 = mat2.ncols();
    let values1 = column_major_copy(mat1.data(), n, p1, center);
    let values2 = column_major_copy(mat2.data(), n, p2, center);
    let denom = n as f64 - 1.0;

    let mut out = RMatrix::<f64>::new(p1, p2);
    let cov = out.data_mut();
    crossprod_blocks(&values1, &values2, n, p1, p2, cov, |_, _| p1);
    cov.par_iter_mut().for_each(|v| *v /= denom);
    Ok(out)
}

pub fn fast_rbind_impl(
    mat1: &RMatrix<f64>,
    mat2: &RMatrix<f64>,
) -> extendr_api::Result<RMatrix<f64>> {
    let ncols = mat1.ncols();
    if mat2.ncols() != ncols {
        return Err(extendr_api::Error::Other(format!(
            "FastRBind: mat1 has {ncols} columns but mat2 has {}.",
            mat2.ncols()
        )));
    }
    let n1 = mat1.nrows();
    let n2 = mat2.nrows();
    let n = n1 + n2;
    let a = mat1.data();
    let b = mat2.data();

    let mut out = RMatrix::<f64>::new(n, ncols);
    if n > 0 {
        for (c, col) in out.data_mut().chunks_mut(n).enumerate() {
            col[..n1].copy_from_slice(&a[c * n1..(c + 1) * n1]);
            col[n1..].copy_from_slice(&b[c * n2..(c + 1) * n2]);
        }
    }
    Ok(out)
}

pub fn row_var_impl(mat: &RMatrix<f64>) -> Doubles {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    let values = mat.data();
    let denom = ncols as f64 - 1.0;

    let out: Vec<f64> = (0..nrows)
        .into_par_iter()
        .map(|r| {
            let row_mean = if ncols == 0 {
                0.0
            } else {
                (0..ncols).map(|c| values[r + c * nrows]).sum::<f64>() / ncols as f64
            };
            (0..ncols)
                .map(|c| {
                    let d = values[r + c * nrows] - row_mean;
                    d * d
                })
                .sum::<f64>()
                / denom
        })
        .collect();

    Doubles::from_values(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_crossprod(a: &[f64], b: &[f64], n: usize, m: usize, q: usize) -> Vec<f64> {
        let mut out = vec![0.0; m * q];
        for j in 0..q {
            for i in 0..m {
                out[i + j * m] = (0..n).map(|k| a[k + i * n] * b[k + j * n]).sum();
            }
        }
        out
    }

    fn patterned(len: usize, seed: u64) -> Vec<f64> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5
            })
            .collect()
    }

    #[test]
    fn crossprod_blocks_matches_naive_product() {
        let (n, m, q) = (37, 70, 101);
        let a = patterned(n * m, 3);
        let b = patterned(n * q, 5);
        let mut out = vec![0.0; m * q];
        crossprod_blocks(&a, &b, n, m, q, &mut out, |_, _| m);
        let expected = naive_crossprod(&a, &b, n, m, q);
        for (g, e) in out.iter().zip(expected.iter()) {
            assert!((g - e).abs() < 1e-12, "{g} vs {e}");
        }
    }

    #[test]
    fn crossprod_blocks_upper_triangle_is_exact_where_computed() {
        let (n, p) = (25, 90);
        let a = patterned(n * p, 11);
        let mut out = vec![0.0; p * p];
        crossprod_blocks(&a, &a, n, p, p, &mut out, |_, j1| j1);
        let expected = naive_crossprod(&a, &a, n, p, p);
        for j in 0..p {
            for i in 0..=j {
                assert!((out[i + j * p] - expected[i + j * p]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn column_major_copy_centers_columns() {
        let values = vec![1.0, 2.0, 3.0, 10.0, 20.0, 30.0];
        let centered = column_major_copy(&values, 3, 2, true);
        assert_eq!(centered, vec![-1.0, 0.0, 1.0, -10.0, 0.0, 10.0]);
        assert_eq!(column_major_copy(&values, 3, 2, false), values);
    }
}
