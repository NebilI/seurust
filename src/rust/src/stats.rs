use extendr_api::prelude::*;

pub fn row_sum_dgcmatrix_impl(x: Doubles, i: Integers, rows: i32) -> Doubles {
    let n_rows = rows as usize;
    let mut rowsum = vec![0.0_f64; n_rows];
    for k in 0..x.len() {
        rowsum[i[k].0 as usize] += x[k].0;
    }
    Doubles::from_values(rowsum)
}

pub fn row_mean_dgcmatrix_impl(x: Doubles, i: Integers, rows: i32, cols: i32) -> Doubles {
    let rowsum = row_sum_dgcmatrix_impl(x, i, rows);
    let n_cols = cols as f64;
    let values: Vec<f64> = (0..rowsum.len()).map(|k| rowsum[k].0 / n_cols).collect();
    Doubles::from_values(values)
}

pub fn row_var_dgcmatrix_impl(x: Doubles, i: Integers, rows: i32, cols: i32) -> Doubles {
    let rowmean = row_mean_dgcmatrix_impl(x.clone(), i.clone(), rows, cols);
    let n_rows = rows as usize;
    let n_cols = cols as i32;

    let mut rowvar = vec![0.0_f64; n_rows];
    let mut nzero = vec![n_cols; n_rows];

    for k in 0..x.len() {
        let row_idx = i[k].0 as usize;
        let diff = x[k].0 - rowmean[row_idx].0;
        rowvar[row_idx] += diff * diff;
        nzero[row_idx] -= 1;
    }

    let denom = (n_cols - 1) as f64;
    for k in 0..n_rows {
        rowvar[k] =
            (rowvar[k] + (rowmean[k].0 * rowmean[k].0 * nzero[k] as f64)) / denom;
    }

    Doubles::from_values(rowvar)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_sum_matches_manual() {
        let x = Doubles::from_values(vec![1.0, 2.0, 3.0, 4.0]);
        let i = Integers::from_values(vec![0, 1, 0, 2]);
        let result = row_sum_dgcmatrix_impl(x, i, 3);
        assert_eq!(result.iter().map(|v| v.0).collect::<Vec<_>>(), vec![4.0, 2.0, 4.0]);
    }

    #[test]
    fn row_mean_divides_by_cols() {
        let x = Doubles::from_values(vec![2.0, 4.0]);
        let i = Integers::from_values(vec![0, 1]);
        let result = row_mean_dgcmatrix_impl(x, i, 2, 4);
        assert_eq!(result.iter().map(|v| v.0).collect::<Vec<_>>(), vec![0.5, 1.0]);
    }
}
