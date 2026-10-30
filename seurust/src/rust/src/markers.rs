//! FindAllMarkers sparse fold-change stats (Seurat 5.6).

use crate::sparse::CscView;
use extendr_api::prelude::*;
use rayon::prelude::*;

fn dimmed_matrix(values: Vec<f64>, nrow: usize, ncol: usize) -> RMatrix<f64> {
    crate::sparse::rmatrix_from_column_major(&values, nrow, ncol)
}

fn dimmed_integers(values: Vec<i32>, nrow: i32, ncol: i32) -> Robj {
    let mut obj = Robj::from(Integers::from_values(values));
    let _ = obj.set_attrib("dim", Integers::from_values(vec![nrow, ncol]));
    obj
}

pub fn find_all_markers_sparse_fold_change_stats_impl(
    view: CscView<'_>,
    groups: &[i32],
    n_groups: i32,
    log_normalize: bool,
) -> extendr_api::Result<Robj> {
    let rows = view.nrows as usize;
    let cols = view.ncols as usize;
    if groups.len() != cols {
        return Err(extendr_api::Error::Other(
            "Length of groups must match the number of matrix columns.".into(),
        ));
    }
    let n_groups = n_groups.max(0) as usize;
    let matrix_size = rows * n_groups;

    let (group_sums, detected) = (0..cols)
        .into_par_iter()
        .fold(
            || (vec![0.0; matrix_size], vec![0i32; matrix_size]),
            |mut acc, col| {
                let group = groups[col] - 1;
                if group < 0 || group as usize >= n_groups {
                    return acc;
                }
                let group_offset = group as usize * rows;
                for ptr in view.p[col] as usize..view.p[col + 1] as usize {
                    let row = view.i[ptr] as usize;
                    let value = view.x[ptr];
                    let sum_value = if log_normalize { value.exp_m1() } else { value };
                    acc.0[group_offset + row] += sum_value;
                    if value > 0.0 {
                        acc.1[group_offset + row] += 1;
                    }
                }
                acc
            },
        )
        .reduce(
            || (vec![0.0; matrix_size], vec![0i32; matrix_size]),
            |mut left, right| {
                for idx in 0..matrix_size {
                    left.0[idx] += right.0[idx];
                    left.1[idx] += right.1[idx];
                }
                left
            },
        );

    let mut rest_sums = vec![0.0; matrix_size];
    let mut total_detected = vec![0i32; rows];
    for row in 0..rows {
        let mut detected_sum = 0i32;
        for group in 0..n_groups {
            detected_sum += detected[row + group * rows];
        }
        total_detected[row] = detected_sum;
        for group in 0..n_groups {
            // Kahan sum of the other groups, matching the C++ rest_sum path.
            let mut rest = 0.0;
            let mut correction = 0.0;
            for other in 0..n_groups {
                if other == group {
                    continue;
                }
                let y = group_sums[row + other * rows] - correction;
                let t = rest + y;
                correction = (t - rest) - y;
                rest = t;
            }
            rest_sums[row + group * rows] = rest;
        }
    }

    Ok(Robj::from(list!(
        group_sum = dimmed_matrix(group_sums, rows, n_groups),
        rest_sum = dimmed_matrix(rest_sums, rows, n_groups),
        detected = dimmed_integers(detected, rows as i32, n_groups as i32),
        total_detected = Integers::from_values(total_detected)
    )))
}
