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

pub fn fast_dist_impl(x: &RMatrix<f64>, y: &RMatrix<f64>, n: &List) -> Robj {
    let ngraph_size = n.len();
    if x.nrows() != ngraph_size {
        return Robj::from(List::new(0));
    }

    let ncols = x.ncols();
    let nrows_x = x.nrows();
    let nrows_y = y.nrows();
    let x_data = x.as_robj().as_real_slice().expect("numeric x");
    let y_data = y.as_robj().as_real_slice().expect("numeric y");

    let neighbors_by_row: Vec<Vec<usize>> = (0..ngraph_size)
        .map(|i| {
            let neighbors: Doubles = n.elt(i).unwrap().try_into().unwrap();
            neighbors.iter().map(|idx| idx.0 as usize - 1).collect()
        })
        .collect();

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

    let mut items = Vec::with_capacity(ngraph_size);
    for distances in distances_by_row {
        items.push(Robj::from(Doubles::from_values(distances)));
    }

    Robj::from(items)
}
