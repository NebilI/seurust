use rayon::prelude::*;

/// Primitive-op count below which rayon dispatch usually costs more than the work.
pub const SERIAL_WORK_THRESHOLD: usize = 65_536;

#[inline]
pub fn use_parallel(n_tasks: usize, work: usize) -> bool {
    n_tasks > 1 && work >= SERIAL_WORK_THRESHOLD
}

/// Column ranges `[bounds[k], bounds[k + 1])` that split `ncols` into
/// `n_chunks` contiguous pieces, balanced by stored nonzeros in CSC `p`.
pub fn col_chunk_bounds(p: &[i32], ncols: usize, n_chunks: usize) -> Vec<usize> {
    let n_chunks = n_chunks.max(1).min(ncols.max(1));
    let mut bounds = vec![0usize; n_chunks + 1];
    bounds[n_chunks] = ncols;
    if n_chunks == 1 || ncols == 0 {
        return bounds;
    }
    let nnz = p.get(ncols).copied().unwrap_or(0).max(0) as usize;
    let mut col = 0;
    for k in 1..n_chunks {
        let target = ((nnz * k) / n_chunks) as i32;
        while col < ncols && p[col] < target {
            col += 1;
        }
        bounds[k] = col;
    }
    bounds
}

/// Run `work(col_start, col_end)` once per nnz-balanced column chunk.
///
/// One chunk per rayon thread (or a single serial call when the matrix is
/// small), so each accumulator is allocated once per thread rather than once
/// per work-stealing split.
pub fn par_col_chunk_map<T, F>(p: &[i32], ncols: usize, work: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize, usize) -> T + Sync + Send,
{
    if ncols == 0 {
        return Vec::new();
    }
    let nnz = p.get(ncols).copied().unwrap_or(0).max(0) as usize;
    let n_chunks = if use_parallel(ncols, nnz) {
        rayon::current_num_threads().max(1).min(ncols)
    } else {
        1
    };
    let bounds = col_chunk_bounds(p, ncols, n_chunks);
    let n = bounds.len() - 1;
    if n <= 1 {
        vec![work(0, ncols)]
    } else {
        (0..n)
            .into_par_iter()
            .map(|k| work(bounds[k], bounds[k + 1]))
            .collect()
    }
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

    #[test]
    fn col_chunk_bounds_cover_every_column_once() {
        let p = vec![0, 3, 5, 5, 10, 12];
        let bounds = col_chunk_bounds(&p, 5, 3);
        assert_eq!(bounds[0], 0);
        assert_eq!(*bounds.last().unwrap(), 5);
        for w in bounds.windows(2) {
            assert!(w[0] <= w[1]);
        }
        let serial = col_chunk_bounds(&p, 5, 1);
        assert_eq!(serial, vec![0, 5]);
        assert_eq!(col_chunk_bounds(&p, 0, 4), vec![0, 0]);
    }
}
