//! SCTransform Pearson residuals and corrected counts (Seurat 5.6).

use crate::sparse::{rmatrix_with_data, CscSlots, CscView};
use extendr_api::prelude::*;
use rayon::prelude::*;

fn sct_model_var(mu: f64, theta: f64) -> f64 {
    if theta.is_finite() {
        mu + (mu * mu / theta)
    } else {
        mu
    }
}

fn sct_clip(value: f64, lower: f64, upper: f64) -> f64 {
    value.max(lower).min(upper)
}

fn residual_stats(
    view: CscView<'_>,
    theta: &[f64],
    intercept: &[f64],
    slope: &[f64],
    log_umi: &[f64],
    min_var: f64,
    residual_clip_min: f64,
    residual_clip_max: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, bool, f64) {
    let rows = view.nrows as usize;
    let cols = view.ncols as usize;
    let exp_intercept: Vec<f64> = intercept.iter().map(|v| v.exp()).collect();
    let first_slope = if rows == 0 { 0.0 } else { slope[0] };
    let common_slope = slope.iter().all(|s| (s - first_slope).abs() <= 1e-12);

    let (residual_sum, residual_sq_sum) = (0..cols)
        .into_par_iter()
        .fold(
            || (vec![0.0; rows], vec![0.0; rows]),
            |mut acc, col| {
                let log_umi_col = log_umi[col];
                let common_factor = if common_slope {
                    (first_slope * log_umi_col).exp()
                } else {
                    0.0
                };
                let start = view.p[col] as usize;
                let end = view.p[col + 1] as usize;
                let mut ptr = start;
                for row in 0..rows {
                    let mut y_value = 0.0;
                    if ptr < end && view.i[ptr] as usize == row {
                        y_value = view.x[ptr];
                        ptr += 1;
                    }
                    let mu_original = if common_slope {
                        exp_intercept[row] * common_factor
                    } else {
                        exp_intercept[row] * (slope[row] * log_umi_col).exp()
                    };
                    let mut variance_original = sct_model_var(mu_original, theta[row]);
                    if variance_original < min_var {
                        variance_original = min_var;
                    }
                    let residual = (y_value - mu_original) / variance_original.sqrt();
                    let clipped = sct_clip(residual, residual_clip_min, residual_clip_max);
                    acc.0[row] += clipped;
                    acc.1[row] += clipped * clipped;
                }
                acc
            },
        )
        .reduce(
            || (vec![0.0; rows], vec![0.0; rows]),
            |mut left, right| {
                for row in 0..rows {
                    left.0[row] += right.0[row];
                    left.1[row] += right.1[row];
                }
                left
            },
        );

    (
        residual_sum,
        residual_sq_sum,
        exp_intercept,
        slope.to_vec(),
        common_slope,
        first_slope,
    )
}

fn corrected_counts(
    view: CscView<'_>,
    theta: &[f64],
    slope: &[f64],
    log_umi: &[f64],
    exp_intercept: &[f64],
    target_log_umi: f64,
    min_var: f64,
    common_slope: bool,
    first_slope: f64,
) -> CscSlots {
    let rows = view.nrows as usize;
    let cols = view.ncols as usize;
    let mut target_mu = vec![0.0; rows];
    let mut target_sqrt_var = vec![0.0; rows];
    if common_slope {
        let target_factor = (first_slope * target_log_umi).exp();
        for row in 0..rows {
            target_mu[row] = exp_intercept[row] * target_factor;
            target_sqrt_var[row] = sct_model_var(target_mu[row], theta[row]).sqrt();
        }
    } else {
        for row in 0..rows {
            target_mu[row] = exp_intercept[row] * (slope[row] * target_log_umi).exp();
            target_sqrt_var[row] = sct_model_var(target_mu[row], theta[row]).sqrt();
        }
    }

    let mut corrected_i = Vec::with_capacity(view.x.len());
    let mut corrected_x = Vec::with_capacity(view.x.len());
    let mut corrected_p = vec![0i32; cols + 1];
    for col in 0..cols {
        corrected_p[col] = corrected_i.len() as i32;
        let start = view.p[col] as usize;
        let end = view.p[col + 1] as usize;
        let mut ptr = start;
        let log_umi_col = log_umi[col];
        let common_factor = if common_slope {
            (first_slope * log_umi_col).exp()
        } else {
            0.0
        };
        for row in 0..rows {
            let mut y_value = 0.0;
            if ptr < end && view.i[ptr] as usize == row {
                y_value = view.x[ptr];
                ptr += 1;
            }
            let mu_original = if common_slope {
                exp_intercept[row] * common_factor
            } else {
                exp_intercept[row] * (slope[row] * log_umi_col).exp()
            };
            let mut variance_original = sct_model_var(mu_original, theta[row]);
            if variance_original < min_var {
                variance_original = min_var;
            }
            let residual = (y_value - mu_original) / variance_original.sqrt();
            let mut corrected = target_mu[row] + residual * target_sqrt_var[row];
            corrected = corrected.round_ties_even();
            if corrected < 0.0 {
                corrected = 0.0;
            }
            if corrected != 0.0 {
                corrected_i.push(row as i32);
                corrected_x.push(corrected);
            }
        }
    }
    corrected_p[cols] = corrected_i.len() as i32;
    CscSlots {
        x: corrected_x,
        i: corrected_i,
        p: corrected_p,
        nrows: view.nrows,
        ncols: view.ncols,
    }
}

pub fn sct_residual_stats_and_corrected_impl(
    view: CscView<'_>,
    theta: &[f64],
    intercept: &[f64],
    slope: &[f64],
    log_umi: &[f64],
    target_log_umi: f64,
    min_var: f64,
    residual_clip_min: f64,
    residual_clip_max: f64,
    compute_corrected: bool,
) -> extendr_api::Result<Robj> {
    let rows = view.nrows as usize;
    let cols = view.ncols as usize;
    let (residual_sum, residual_sq_sum, exp_intercept, slope_owned, common_slope, first_slope) =
        residual_stats(
            view,
            theta,
            intercept,
            slope,
            log_umi,
            min_var,
            residual_clip_min,
            residual_clip_max,
        );
    let cols_d = cols as f64;
    let denom = (cols.max(2) - 1) as f64;
    let residual_mean: Vec<f64> = residual_sum.iter().map(|s| s / cols_d).collect();
    let residual_variance: Vec<f64> = (0..rows)
        .map(|row| (residual_sq_sum[row] - residual_sum[row] * residual_sum[row] / cols_d) / denom)
        .collect();

    if compute_corrected {
        let corrected = corrected_counts(
            view,
            theta,
            &slope_owned,
            log_umi,
            &exp_intercept,
            target_log_umi,
            min_var,
            common_slope,
            first_slope,
        );
        Ok(Robj::from(list!(
            residual_mean = Doubles::from_values(residual_mean),
            residual_variance = Doubles::from_values(residual_variance),
            corrected = corrected.into_r_dgcmatrix()?
        )))
    } else {
        Ok(Robj::from(list!(
            residual_mean = Doubles::from_values(residual_mean),
            residual_variance = Doubles::from_values(residual_variance)
        )))
    }
}

pub fn sct_pearson_residual_matrix_impl(
    view: CscView<'_>,
    theta: &[f64],
    intercept: &[f64],
    slope: &[f64],
    log_umi: &[f64],
    feature_index: &[i32],
    min_var: &[f64],
    clip_min: f64,
    clip_max: f64,
    do_center: bool,
) -> extendr_api::Result<RMatrix<f64>> {
    let rows = view.nrows as usize;
    let cols = view.ncols as usize;
    let selected = feature_index.len();
    if min_var.len() != 1 && min_var.len() != selected {
        return Err(extendr_api::Error::Other(
            "min_var must have length 1 or match feature_index".into(),
        ));
    }
    let mut row_to_selected = vec![-1i32; rows];
    let mut selected_theta = vec![0.0; selected];
    let mut selected_exp_intercept = vec![0.0; selected];
    let mut selected_slope = vec![0.0; selected];
    let mut selected_min_var = vec![0.0; selected];
    let first_slope = if selected > 0 {
        slope[feature_index[0] as usize]
    } else {
        0.0
    };
    let mut common_slope = true;
    for (idx, &feat) in feature_index.iter().enumerate() {
        let row = feat as usize;
        row_to_selected[row] = idx as i32;
        selected_theta[idx] = theta[row];
        selected_exp_intercept[idx] = intercept[row].exp();
        selected_slope[idx] = slope[row];
        selected_min_var[idx] = if min_var.len() == 1 {
            min_var[0]
        } else {
            min_var[idx]
        };
        if (selected_slope[idx] - first_slope).abs() > 1e-12 {
            common_slope = false;
        }
    }

    let mut row_sum = vec![0.0; selected];
    let mut out = rmatrix_with_data(selected, cols, |out| {
        let out_addr = out.as_mut_ptr() as usize;
        let sums: Vec<Vec<f64>> = (0..cols)
            .into_par_iter()
            .map(|col| {
                let mut local_y = vec![0.0; selected];
                for ptr in view.p[col] as usize..view.p[col + 1] as usize {
                    let selected_row = row_to_selected[view.i[ptr] as usize];
                    if selected_row >= 0 {
                        local_y[selected_row as usize] = view.x[ptr];
                    }
                }
                let log_umi_col = log_umi[col];
                let common_factor = if common_slope {
                    (first_slope * log_umi_col).exp()
                } else {
                    0.0
                };
                let mut local_sum = vec![0.0; selected];
                let out_ptr = out_addr as *mut f64;
                unsafe {
                    let out_col = out_ptr.add(col * selected);
                    for idx in 0..selected {
                        let mu = if common_slope {
                            selected_exp_intercept[idx] * common_factor
                        } else {
                            selected_exp_intercept[idx] * (selected_slope[idx] * log_umi_col).exp()
                        };
                        let mut variance = sct_model_var(mu, selected_theta[idx]);
                        if variance < selected_min_var[idx] {
                            variance = selected_min_var[idx];
                        }
                        let residual = (local_y[idx] - mu) / variance.sqrt();
                        let clipped = sct_clip(residual, clip_min, clip_max);
                        *out_col.add(idx) = clipped;
                        local_sum[idx] = clipped;
                    }
                }
                local_sum
            })
            .collect();
        for local in sums {
            for idx in 0..selected {
                row_sum[idx] += local[idx];
            }
        }
    });

    if do_center && cols > 0 {
        let inv_cols = 1.0 / cols as f64;
        out.data_mut()
            .par_chunks_mut(selected.max(1))
            .for_each(|col| {
                for idx in 0..selected {
                    col[idx] -= row_sum[idx] * inv_cols;
                }
            });
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn model_var_poisson_when_theta_is_inf() {
        assert_eq!(sct_model_var(4.0, f64::INFINITY), 4.0);
        assert!((sct_model_var(4.0, 2.0) - (4.0 + 8.0)).abs() < 1e-12);
    }
}
