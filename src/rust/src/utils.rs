/// Stable sort of indices by ascending values (matches C++ `sort_indexes`).
pub fn sort_indexes(v: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    idx
}

pub fn euclidean_rows(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

pub fn col_euclidean_dist(mat: &ndarray::Array2<f64>, col_a: usize, col_b: usize) -> f64 {
    let nrows = mat.nrows();
    let mut sum = 0.0;
    for r in 0..nrows {
        let d = mat[[r, col_a]] - mat[[r, col_b]];
        sum += d * d;
    }
    sum.sqrt()
}

pub fn row_euclidean_dist(mat: &ndarray::Array2<f64>, row_a: usize, row_b: usize) -> f64 {
    let ncols = mat.ncols();
    let mut sum = 0.0;
    for c in 0..ncols {
        let d = mat[[row_a, c]] - mat[[row_b, c]];
        sum += d * d;
    }
    sum.sqrt()
}

pub fn csc_to_dense(slots: &crate::sparse::CscSlots) -> ndarray::Array2<f64> {
    let nrows = slots.nrows as usize;
    let ncols = slots.ncols as usize;
    let mut out = ndarray::Array2::zeros((nrows, ncols));
    for col in 0..ncols {
        for idx in slots.p[col] as usize..slots.p[col + 1] as usize {
            out[[slots.i[idx] as usize, col]] = slots.x[idx];
        }
    }
    out
}

pub fn dense_to_csc(mat: &ndarray::Array2<f64>) -> crate::sparse::CscSlots {
    let (nrows, ncols) = mat.dim();
    let mut triplets = Vec::new();
    for r in 0..nrows {
        for c in 0..ncols {
            let v = mat[[r, c]];
            if v != 0.0 {
                triplets.push((r, c, v));
            }
        }
    }
    crate::sparse::csc_from_triplets(nrows, ncols, &triplets)
}
