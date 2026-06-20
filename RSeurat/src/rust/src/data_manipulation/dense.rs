use crate::sparse::{ndarray_from_rmatrix, rmatrix_from_column_major, rmatrix_from_ndarray};
use extendr_api::prelude::*;
use ndarray::Axis;
use rayon::prelude::*;

pub fn standardize_impl(mat: &RMatrix<f64>, _display_progress: bool) -> RMatrix<f64> {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    let values = mat.as_robj().as_real_slice().expect("numeric mat");
    let mut out = vec![0.0; nrows * ncols];

    out.par_chunks_mut(nrows)
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
                / (nrows - 1) as f64)
                .sqrt();
            for r in 0..nrows {
                out_col[r] = (col[r] - col_mean) / col_sdev;
            }
        });

    rmatrix_from_column_major(&out, nrows, ncols)
}

pub fn fast_cov_impl(mat: &RMatrix<f64>, center: bool) -> RMatrix<f64> {
    let mut values = ndarray_from_rmatrix(mat);
    if center {
        let means = values.mean_axis(Axis(0)).unwrap();
        for mut row in values.rows_mut() {
            row -= &means;
        }
    }
    let nrows = values.nrows() as f64;
    let cov = values.t().dot(&values) / (nrows - 1.0);
    rmatrix_from_ndarray(cov.view())
}

pub fn fast_cov_mats_impl(mat1: &RMatrix<f64>, mat2: &RMatrix<f64>, center: bool) -> RMatrix<f64> {
    let mut values1 = ndarray_from_rmatrix(mat1);
    let mut values2 = ndarray_from_rmatrix(mat2);
    if center {
        let means1 = values1.mean_axis(Axis(0)).unwrap();
        let means2 = values2.mean_axis(Axis(0)).unwrap();
        for mut row in values1.rows_mut() {
            row -= &means1;
        }
        for mut row in values2.rows_mut() {
            row -= &means2;
        }
    }
    let nrows = values1.nrows() as f64;
    let cov = values1.t().dot(&values2) / (nrows - 1.0);
    rmatrix_from_ndarray(cov.view())
}

pub fn fast_rbind_impl(mat1: &RMatrix<f64>, mat2: &RMatrix<f64>) -> RMatrix<f64> {
    let values1 = ndarray_from_rmatrix(mat1);
    let values2 = ndarray_from_rmatrix(mat2);
    let combined = ndarray::concatenate![Axis(0), values1, values2];
    rmatrix_from_ndarray(combined.view())
}

pub fn row_var_impl(mat: &RMatrix<f64>) -> Doubles {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    let values = mat.as_robj().as_real_slice().expect("numeric mat");
    let denom = (ncols - 1) as f64;

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
