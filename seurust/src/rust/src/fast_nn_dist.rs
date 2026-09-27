use extendr_api::prelude::*;
use rayon::prelude::*;

fn transpose_row_major(data: &[f64], nrows: usize, ncols: usize) -> Vec<f64> {
    let mut out = vec![0.0; nrows.saturating_mul(ncols)];
    for c in 0..ncols {
        let col = &data[c * nrows..(c + 1) * nrows];
        for r in 0..nrows {
            out[r * ncols + c] = col[r];
        }
    }
    out
}

fn euclidean(xs: &[f64], ys: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (&a, &b) in xs.iter().zip(ys) {
        let d = a - b;
        sum += d * d;
    }
    sum.sqrt()
}

/// Convert one element of the neighbor list to 0-based row indices into `y`.
/// Seurat's Rcpp signature coerces integer and double vectors alike.
fn neighbor_indices(elt: &Robj, cell: usize, nrows_y: usize) -> extendr_api::Result<Vec<usize>> {
    let to_index = |raw: f64| -> extendr_api::Result<usize> {
        if !raw.is_finite() || raw < 1.0 || raw > nrows_y as f64 {
            return Err(extendr_api::Error::Other(format!(
                "fast_dist: neighbor index {raw} for cell {} is outside 1..{nrows_y}.",
                cell + 1
            )));
        }
        Ok(raw as usize - 1)
    };

    if let Some(ints) = elt.as_integer_slice() {
        ints.iter()
            .map(|&v| {
                if v == i32::MIN {
                    to_index(f64::NAN)
                } else {
                    to_index(v as f64)
                }
            })
            .collect()
    } else if let Some(reals) = elt.as_real_slice() {
        reals.iter().map(|&v| to_index(v)).collect()
    } else {
        Err(extendr_api::Error::Other(format!(
            "fast_dist: neighbors for cell {} must be an integer or numeric vector.",
            cell + 1
        )))
    }
}

fn neighbor_id(raw: f64, cell: usize, nrows_y: usize) -> extendr_api::Result<usize> {
    if !raw.is_finite() || raw < 1.0 || raw > nrows_y as f64 {
        return Err(extendr_api::Error::Other(format!(
            "fast_dist: neighbor index {raw} for cell {} is outside 1..{nrows_y}.",
            cell + 1
        )));
    }
    Ok(raw as usize - 1)
}

fn write_distances_column_major(
    out: &mut [f64],
    elt: &Robj,
    x: &[f64],
    y: &[f64],
    x_row: usize,
    ncols: usize,
    nrows_x: usize,
    nrows_y: usize,
) -> extendr_api::Result<()> {
    let fill = |out: &mut [f64], n_idx: usize, j: usize| {
        let mut sum = 0.0;
        for c in 0..ncols {
            let d = x[x_row + c * nrows_x] - y[n_idx + c * nrows_y];
            sum += d * d;
        }
        out[j] = sum.sqrt();
    };
    let cell = x_row;

    if let Some(ints) = elt.as_integer_slice() {
        if ints.len() != out.len() {
            return Err(extendr_api::Error::Other(
                "fast_dist: internal length mismatch".to_string(),
            ));
        }
        for (j, &v) in ints.iter().enumerate() {
            let raw = if v == i32::MIN { f64::NAN } else { v as f64 };
            fill(out, neighbor_id(raw, cell, nrows_y)?, j);
        }
        Ok(())
    } else if let Some(reals) = elt.as_real_slice() {
        if reals.len() != out.len() {
            return Err(extendr_api::Error::Other(
                "fast_dist: internal length mismatch".to_string(),
            ));
        }
        for (j, &raw) in reals.iter().enumerate() {
            fill(out, neighbor_id(raw, cell, nrows_y)?, j);
        }
        Ok(())
    } else {
        Err(extendr_api::Error::Other(format!(
            "fast_dist: neighbors for cell {} must be an integer or numeric vector.",
            cell + 1
        )))
    }
}

pub fn fast_dist_impl(x: &RMatrix<f64>, y: &RMatrix<f64>, n: &List) -> extendr_api::Result<Robj> {
    let ngraph_size = n.len();
    if x.nrows() != ngraph_size {
        return Ok(Robj::from(List::new(0)));
    }
    if x.ncols() != y.ncols() {
        return Err(extendr_api::Error::Other(format!(
            "fast_dist: x has {} columns but y has {}.",
            x.ncols(),
            y.ncols()
        )));
    }

    let ncols = x.ncols();
    let nrows_y = y.nrows();
    let x_data = x.data();
    let y_data = y.data();

    // R object allocation stays on this thread. Parallelism pays off once the
    // distance work dwarfs scheduling; the published 800 x 20 case does not.
    let items: Vec<Robj> = if ngraph_size >= 4096 {
        let x_rm = transpose_row_major(x_data, ngraph_size, ncols);
        let y_rm = transpose_row_major(y_data, nrows_y, ncols);
        let neighbors_by_row = n
            .values()
            .enumerate()
            .map(|(i, elt)| neighbor_indices(&elt, i, nrows_y))
            .collect::<extendr_api::Result<Vec<Vec<usize>>>>()?;
        let distances_by_row: Vec<Vec<f64>> = neighbors_by_row
            .par_iter()
            .enumerate()
            .map(|(i, neighbors)| {
                let xs = &x_rm[i * ncols..(i + 1) * ncols];
                neighbors
                    .iter()
                    .map(|&n_idx| euclidean(xs, &y_rm[n_idx * ncols..(n_idx + 1) * ncols]))
                    .collect()
            })
            .collect();
        distances_by_row
            .into_iter()
            .map(|distances| Robj::from(Doubles::from_values(distances)))
            .collect()
    } else {
        let mut items = Vec::with_capacity(ngraph_size);
        for (i, elt) in n.values().enumerate() {
            let len = if let Some(ints) = elt.as_integer_slice() {
                ints.len()
            } else if let Some(reals) = elt.as_real_slice() {
                reals.len()
            } else {
                return Err(extendr_api::Error::Other(format!(
                    "fast_dist: neighbors for cell {} must be an integer or numeric vector.",
                    i + 1
                )));
            };
            let mut distances = Doubles::new(len);
            if len > 0 {
                write_distances_column_major(
                    distances
                        .as_robj_mut()
                        .as_real_slice_mut()
                        .expect("numeric distances"),
                    &elt,
                    x_data,
                    y_data,
                    i,
                    ncols,
                    ngraph_size,
                    nrows_y,
                )?;
            }
            items.push(Robj::from(distances));
        }
        items
    };

    let mut out = Robj::from(List::from_values(items));
    if let Some(names) = n.names() {
        let names: Vec<&str> = names.collect();
        out.set_names(names)?;
    }
    Ok(out)
}
