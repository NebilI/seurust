use crate::row_sum_kernels::row_sums_into;
use crate::sparse::doubles_uninit;
use extendr_api::prelude::*;

pub fn row_sum_slices(x: &[f64], i: &[i32], rows: i32) -> Vec<f64> {
    let n_rows = rows as usize;
    let mut rowsum = vec![0.0; n_rows];
    for (&row, &val) in i.iter().zip(x.iter()) {
        rowsum[row as usize] += val;
    }
    rowsum
}

/// Row sums (or means when `divisor` is set) written straight into a fresh R
/// vector: no intermediate `Vec`, no copy, no separate divide pass. Large
/// inputs are split across rayon threads (see `row_sum_kernels`).
fn row_sums_to_r(x: &Doubles, i: &Integers, rows: i32, divisor: Option<f64>) -> Doubles {
    let n_rows = usize::try_from(rows).expect("non-negative rows");
    let x_data = x.as_robj().as_real_slice().expect("numeric x");
    let i_data = i.as_robj().as_integer_slice().expect("integer i");
    let mut out = doubles_uninit(n_rows);
    if n_rows > 0 {
        let buf = out
            .as_robj_mut()
            .as_real_slice_mut()
            .expect("numeric output");
        row_sums_into(x_data, i_data, buf, divisor);
    }
    out
}

pub fn row_sum_dgcmatrix_impl(x: &Doubles, i: &Integers, rows: i32) -> Doubles {
    row_sums_to_r(x, i, rows, None)
}

pub fn row_mean_dgcmatrix_impl(x: &Doubles, i: &Integers, rows: i32, cols: i32) -> Doubles {
    row_sums_to_r(x, i, rows, Some(cols as f64))
}

pub fn row_var_dgcmatrix_impl(x: &Doubles, i: &Integers, rows: i32, cols: i32) -> Doubles {
    let n_rows = rows as usize;
    let n_cols = cols as i32;
    let ncol_f = cols as f64;
    let denom = (n_cols - 1) as f64;

    let x_data = x.as_robj().as_real_slice().expect("numeric x");
    let i_data = i.as_robj().as_integer_slice().expect("integer i");

    let rowsum = row_sum_slices(x_data, i_data, rows);
    let mut rowvar = vec![0.0; n_rows];
    let mut nnz = vec![0i32; n_rows];

    for (&row, &val) in i_data.iter().zip(x_data.iter()) {
        let row_idx = row as usize;
        let mean = rowsum[row_idx] / ncol_f;
        let diff = val - mean;
        rowvar[row_idx] += diff * diff;
        nnz[row_idx] += 1;
    }

    for k in 0..n_rows {
        let mean = rowsum[k] / ncol_f;
        let nzero = n_cols - nnz[k];
        rowvar[k] = (rowvar[k] + mean * mean * nzero as f64) / denom;
    }

    Doubles::from_values(rowvar)
}

#[cfg(test)]
mod tests {
    use super::row_sum_slices;

    #[test]
    fn row_sum_slices_adds_by_row_index() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let i = [0, 0, 1, 2];
        let sums = row_sum_slices(&x, &i, 3);
        assert_eq!(sums, vec![3.0, 3.0, 4.0]);
    }

    #[test]
    fn row_sum_slices_handles_empty_input() {
        let sums = row_sum_slices(&[], &[], 2);
        assert_eq!(sums, vec![0.0, 0.0]);
    }

    #[test]
    fn row_sum_matches_manual() {
        let result = row_sum_slices(&[1.0, 2.0, 3.0, 4.0], &[0, 1, 0, 2], 3);
        assert_eq!(result, vec![4.0, 2.0, 4.0]);
    }

    #[test]
    fn row_mean_divides_by_cols() {
        let mut rowsum = row_sum_slices(&[2.0, 4.0], &[0, 1], 2);
        for v in rowsum.iter_mut() {
            *v /= 4.0;
        }
        assert_eq!(rowsum, vec![0.5, 1.0]);
    }
}
