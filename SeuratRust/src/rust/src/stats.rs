use extendr_api::prelude::*;

pub fn row_sum_dgcmatrix_impl(x: &Doubles, i: &Integers, rows: i32) -> Doubles {
    let n_rows = rows as usize;
    let x_data = x.as_robj().as_real_slice().expect("numeric x");
    let i_data = i.as_robj().as_integer_slice().expect("integer i");

    let mut out = Doubles::new(n_rows);
    let rowsum = out.as_robj_mut().as_real_slice_mut().expect("output");

    for idx in 0..x_data.len() {
        rowsum[i_data[idx] as usize] += x_data[idx];
    }
    out
}

pub fn row_mean_dgcmatrix_impl(x: &Doubles, i: &Integers, rows: i32, cols: i32) -> Doubles {
    let mut out = row_sum_dgcmatrix_impl(x, i, rows);
    let n_cols = cols as f64;
    let rowsum = out.as_robj_mut().as_real_slice_mut().expect("output");
    for v in rowsum.iter_mut() {
        *v /= n_cols;
    }
    out
}

pub fn row_var_dgcmatrix_impl(x: &Doubles, i: &Integers, rows: i32, cols: i32) -> Doubles {
    let n_rows = rows as usize;
    let n_cols = cols as i32;
    let ncol_f = cols as f64;
    let denom = (n_cols - 1) as f64;

    let x_data = x.as_robj().as_real_slice().expect("numeric x");
    let i_data = i.as_robj().as_integer_slice().expect("integer i");

    let mut rowsum = vec![0.0_f64; n_rows];
    for (&row, &val) in i_data.iter().zip(x_data.iter()) {
        rowsum[row as usize] += val;
    }

    let mut out = Doubles::new(n_rows);
    let rowvar = out.as_robj_mut().as_real_slice_mut().expect("output");
    rowvar.fill(0.0);

    let mut nzero = vec![n_cols; n_rows];
    for (&row, &val) in i_data.iter().zip(x_data.iter()) {
        let row_idx = row as usize;
        let mean = rowsum[row_idx] / ncol_f;
        let diff = val - mean;
        rowvar[row_idx] += diff * diff;
        nzero[row_idx] -= 1;
    }

    for k in 0..n_rows {
        let mean = rowsum[k] / ncol_f;
        rowvar[k] = (rowvar[k] + mean * mean * nzero[k] as f64) / denom;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_sum_matches_manual() {
        let x = Doubles::from_values(vec![1.0, 2.0, 3.0, 4.0]);
        let i = Integers::from_values(vec![0, 1, 0, 2]);
        let result = row_sum_dgcmatrix_impl(&x, &i, 3);
        assert_eq!(result.iter().map(|v| v.0).collect::<Vec<_>>(), vec![4.0, 2.0, 4.0]);
    }

    #[test]
    fn row_mean_divides_by_cols() {
        let x = Doubles::from_values(vec![2.0, 4.0]);
        let i = Integers::from_values(vec![0, 1]);
        let result = row_mean_dgcmatrix_impl(&x, &i, 2, 4);
        assert_eq!(result.iter().map(|v| v.0).collect::<Vec<_>>(), vec![0.5, 1.0]);
    }
}
