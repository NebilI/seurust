use extendr_api::prelude::*;
use rayon::prelude::*;

fn row_distance_column_major(
    x: &[f64],
    y: &[f64],
    x_row: usize,
    y_row: usize,
    nrows_x: usize,
    nrows_y: usize,
    ncols: usize,
) -> f64 {
    let mut sum = 0.0;
    for c in 0..ncols {
        let d = x[x_row + c * nrows_x] - y[y_row + c * nrows_y];
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
    let nrows_x = x.nrows();
    let nrows_y = y.nrows();
    let x_data = x.data();
    let y_data = y.data();

    let neighbors_by_row = n
        .values()
        .enumerate()
        .map(|(i, elt)| neighbor_indices(&elt, i, nrows_y))
        .collect::<extendr_api::Result<Vec<Vec<usize>>>>()?;

    let distances_by_row: Vec<Vec<f64>> = neighbors_by_row
        .par_iter()
        .enumerate()
        .map(|(i, neighbors)| {
            neighbors
                .iter()
                .map(|&n_idx| {
                    row_distance_column_major(x_data, y_data, i, n_idx, nrows_x, nrows_y, ncols)
                })
                .collect()
        })
        .collect();

    let items: Vec<Robj> = distances_by_row
        .into_iter()
        .map(|distances| Robj::from(Doubles::from_values(distances)))
        .collect();
    let mut out = Robj::from(List::from_values(items));
    if let Some(names) = n.names() {
        let names: Vec<&str> = names.collect();
        out.set_names(names)?;
    }
    Ok(out)
}
