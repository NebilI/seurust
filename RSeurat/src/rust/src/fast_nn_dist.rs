use crate::utils::euclidean_rows;
use crate::sparse::ndarray_from_rmatrix;
use extendr_api::prelude::*;

pub fn fast_dist_impl(x: &RMatrix<f64>, y: &RMatrix<f64>, n: &List) -> Robj {
    let ngraph_size = n.len();
    if x.nrows() != ngraph_size {
        return Robj::from(List::new(0));
    }

    let x_arr = ndarray_from_rmatrix(x);
    let y_arr = ndarray_from_rmatrix(y);
    let mut items = Vec::with_capacity(ngraph_size);

    for i in 0..ngraph_size {
        let neighbors: Doubles = n.elt(i).unwrap().try_into().unwrap();
        let mut distances = Vec::with_capacity(neighbors.len());
        let row_x = x_arr.row(i);

        for j in 0..neighbors.len() {
            let n_idx = neighbors[j].0 as usize - 1;
            let row_y = y_arr.row(n_idx);
            distances.push(euclidean_rows(row_x.as_slice().unwrap(), row_y.as_slice().unwrap()));
        }

        items.push(Robj::from(Doubles::from_values(distances)));
    }

    Robj::from(items)
}
