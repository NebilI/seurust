use crate::sparse::{ndarray_from_rmatrix, rmatrix_from_ndarray};
use extendr_api::prelude::*;
use ndarray::{Array2, Axis};

pub fn standardize_impl(mat: &RMatrix<f64>, _display_progress: bool) -> RMatrix<f64> {
    let values = ndarray_from_rmatrix(mat);
    let (nrows, ncols) = values.dim();
    let mut out = Array2::zeros((nrows, ncols));

    for c in 0..ncols {
        let col = values.column(c);
        let col_mean = col.mean().unwrap_or(0.0);
        let col_sdev = (col.iter().map(|v| (v - col_mean).powi(2)).sum::<f64>() / (nrows - 1) as f64)
            .sqrt();
        for r in 0..nrows {
            out[[r, c]] = (values[[r, c]] - col_mean) / col_sdev;
        }
    }

    rmatrix_from_ndarray(out.view())
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
    let values = ndarray_from_rmatrix(mat);
    let (nrows, ncols) = values.dim();
    let mut out = Vec::with_capacity(nrows);
    let denom = (ncols - 1) as f64;

    for r in 0..nrows {
        let row = values.row(r);
        let row_mean = row.mean().unwrap_or(0.0);
        let var = row.iter().map(|v| (v - row_mean).powi(2)).sum::<f64>() / denom;
        out.push(var);
    }

    Doubles::from_values(out)
}
