use rayon::prelude::*;

/// Primitive-op count below which rayon dispatch usually costs more than the work.
pub const SERIAL_WORK_THRESHOLD: usize = 65_536;

#[inline]
pub fn use_parallel(n_tasks: usize, work: usize) -> bool {
    n_tasks > 1 && work >= SERIAL_WORK_THRESHOLD
}

/// Map `0..n` in parallel only when `work` is large enough to pay for rayon.
pub fn map_indexed<T, F>(n: usize, work: usize, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    if use_parallel(n, work) {
        (0..n).into_par_iter().map(f).collect()
    } else {
        (0..n).map(f).collect()
    }
}

/// Stable sort of indices by ascending values (matches C++ `sort_indexes`).
pub fn sort_indexes(v: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    idx
}

pub fn euclidean_rows(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| {
            let d = x - y;
            d * d
        })
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn map_indexed_serial_matches_parallel() {
        let f = |i: usize| (i as i32).wrapping_mul(3).wrapping_add(1);
        let serial = map_indexed(64, 0, f);
        let parallel = map_indexed(64, SERIAL_WORK_THRESHOLD, f);
        assert_eq!(serial, parallel);
        assert_eq!(serial[7], 22);
        assert!(map_indexed::<i32, _>(0, SERIAL_WORK_THRESHOLD, f).is_empty());
    }

    #[test]
    fn bench_sized_tiny_kernels_stay_serial() {
        assert!(!use_parallel(200, 200 * 80));
        assert!(!use_parallel(80, 80 * 10));
        assert!(!use_parallel(3_200, 80 * 40));
        assert!(use_parallel(2_000, 2_000 * 2_000));
    }
}
