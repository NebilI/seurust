//! Slice-only kernels behind `row_sum_dgcmatrix` and `row_mean_dgcmatrix`.
//!
//! Nothing here touches extendr or R, so the code can be unit-tested and
//! benchmarked without an R installation. The wrappers in `stats.rs` hand in
//! the R output vector as `out`.
//!
//! * [`row_sums_serial`] is the single-pass scatter-add, written straight
//!   into the caller's buffer.
//! * [`row_sums_parallel`] splits the stored entries into equal contiguous
//!   pieces, gives each piece one accumulator, then adds the accumulators up
//!   per row range in parallel.
//! * [`row_sums_into`] picks between them by size.

use rayon::prelude::*;

/// Fewest stored entries a piece must hold before a second thread is worth it.
pub const MIN_NNZ_PER_CHUNK: usize = 100_000;

/// Scatter-add `x` into `out[i]`, then optionally divide every row by `divisor`.
///
/// `out` is overwritten, so it may hold uninitialised R memory. An
/// out-of-range index in `i` panics (extendr turns that into an R error).
pub fn row_sums_serial(x: &[f64], i: &[i32], out: &mut [f64], divisor: Option<f64>) {
    out.fill(0.0);
    for (&row, &val) in i.iter().zip(x) {
        out[row as usize] += val;
    }
    if let Some(d) = divisor {
        for v in out.iter_mut() {
            *v /= d;
        }
    }
}

/// Same result as [`row_sums_serial`] using `n_chunks` threads.
///
/// Chunk 0 accumulates directly into `out`; the others use their own buffer.
/// Each row's total is `((chunk0 + chunk1) + chunk2) ...`. The chunk edges
/// depend only on `x.len()` and `n_chunks`, so the rounding is reproducible
/// for a given thread count. Integer-valued data sums exactly in any order.
pub fn row_sums_parallel(
    x: &[f64],
    i: &[i32],
    out: &mut [f64],
    divisor: Option<f64>,
    n_chunks: usize,
) {
    let nnz = x.len().min(i.len());
    let n_chunks = n_chunks.clamp(1, nnz.max(1));
    if n_chunks == 1 {
        return row_sums_serial(x, i, out, divisor);
    }
    let n_rows = out.len();
    let edge = |k: usize| nnz * k / n_chunks;

    // Phase 1: every chunk scatters its own range of entries.
    let mut extra: Vec<Vec<f64>> = (1..n_chunks).map(|_| vec![0.0; n_rows]).collect();
    out.fill(0.0);
    rayon::scope(|s| {
        let (first, rest) = (&mut *out, &mut extra);
        s.spawn(move |_| {
            for k in edge(0)..edge(1) {
                first[i[k] as usize] += x[k];
            }
        });
        for (c, acc) in rest.iter_mut().enumerate() {
            let (lo, hi) = (edge(c + 1), edge(c + 2));
            s.spawn(move |_| {
                for k in lo..hi {
                    acc[i[k] as usize] += x[k];
                }
            });
        }
    });

    // Phase 2: add the other chunks into `out`, one row range per thread.
    let rows_per_task = n_rows.div_ceil(n_chunks).max(1);
    let extra = &extra;
    out.par_chunks_mut(rows_per_task)
        .enumerate()
        .for_each(|(t, dest)| {
            let r0 = t * rows_per_task;
            for acc in extra {
                for (v, a) in dest.iter_mut().zip(&acc[r0..]) {
                    *v += *a;
                }
            }
            if let Some(d) = divisor {
                for v in dest.iter_mut() {
                    *v /= d;
                }
            }
        });
}

/// Number of pieces to use for `nnz` stored entries on `threads` threads.
pub fn chunk_count(nnz: usize, threads: usize) -> usize {
    (nnz / MIN_NNZ_PER_CHUNK).min(threads).max(1)
}

/// Row sums (or means, when `divisor` is set) into `out`, choosing serial or
/// parallel by the number of stored entries.
pub fn row_sums_into(x: &[f64], i: &[i32], out: &mut [f64], divisor: Option<f64>) {
    let n_chunks = chunk_count(x.len().min(i.len()), rayon::current_num_threads());
    if n_chunks > 1 {
        row_sums_parallel(x, i, out, divisor, n_chunks);
    } else {
        row_sums_serial(x, i, out, divisor);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The previous implementation, kept as the reference.
    fn reference(x: &[f64], i: &[i32], rows: usize, divisor: Option<f64>) -> Vec<f64> {
        let mut s = vec![0.0; rows];
        for (&r, &v) in i.iter().zip(x) {
            s[r as usize] += v;
        }
        if let Some(d) = divisor {
            for v in s.iter_mut() {
                *v /= d;
            }
        }
        s
    }

    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *seed >> 33
    }

    fn input(rows: usize, nnz: usize, integer: bool) -> (Vec<f64>, Vec<i32>) {
        let mut seed = 7;
        let i: Vec<i32> = (0..nnz)
            .map(|_| (lcg(&mut seed) % rows as u64) as i32)
            .collect();
        let x: Vec<f64> = (0..nnz)
            .map(|_| {
                let r = lcg(&mut seed);
                if integer {
                    (1 + r % 20) as f64
                } else {
                    (r % 100_000) as f64 / 997.0
                }
            })
            .collect();
        (x, i)
    }

    #[test]
    fn serial_matches_reference_and_overwrites_garbage() {
        let (x, i) = input(50, 400, false);
        let mut out = vec![f64::NAN; 50];
        row_sums_serial(&x, &i, &mut out, None);
        assert_eq!(out, reference(&x, &i, 50, None));
    }

    #[test]
    fn serial_mean_divides_by_cols() {
        let mut out = vec![9.0; 3];
        row_sums_serial(&[2.0, 4.0, 6.0], &[0, 1, 0], &mut out, Some(4.0));
        assert_eq!(out, vec![2.0, 1.0, 0.0]);
    }

    #[test]
    fn handles_empty_input() {
        let mut out = vec![1.0; 4];
        row_sums_serial(&[], &[], &mut out, None);
        assert_eq!(out, vec![0.0; 4]);
        row_sums_parallel(&[], &[], &mut out, Some(2.0), 3);
        assert_eq!(out, vec![0.0; 4]);
        let mut none: Vec<f64> = vec![];
        row_sums_into(&[], &[], &mut none, None);
    }

    #[test]
    fn parallel_is_bit_identical_on_integer_data() {
        let (x, i) = input(6000, 360_000, true);
        for divisor in [None, Some(1500.0)] {
            let want = reference(&x, &i, 6000, divisor);
            for chunks in 1..=6 {
                let mut out = vec![f64::NAN; 6000];
                row_sums_parallel(&x, &i, &mut out, divisor, chunks);
                assert_eq!(out, want, "chunks = {chunks}, divisor = {divisor:?}");
            }
        }
    }

    #[test]
    fn parallel_is_close_and_deterministic_on_fractional_data() {
        let (x, i) = input(3000, 250_000, false);
        let want = reference(&x, &i, 3000, None);
        let mut a = vec![0.0; 3000];
        let mut b = vec![0.0; 3000];
        row_sums_parallel(&x, &i, &mut a, None, 3);
        row_sums_parallel(&x, &i, &mut b, None, 3);
        assert_eq!(a, b, "same chunk count must give identical bits");
        for (g, w) in a.iter().zip(&want) {
            assert!((g - w).abs() <= 1e-12 * w.abs().max(1.0), "{g} vs {w}");
        }
    }

    #[test]
    fn parallel_handles_more_chunks_than_entries_and_few_rows() {
        let mut out = vec![0.0; 2];
        row_sums_parallel(&[1.0, 2.0, 3.0], &[1, 0, 1], &mut out, None, 16);
        assert_eq!(out, vec![2.0, 4.0]);
    }

    #[test]
    fn dispatch_matches_reference_below_and_above_threshold() {
        for nnz in [1_000, 360_000] {
            let (x, i) = input(6000, nnz, true);
            let mut out = vec![0.0; 6000];
            row_sums_into(&x, &i, &mut out, Some(1500.0));
            assert_eq!(out, reference(&x, &i, 6000, Some(1500.0)));
        }
    }

    #[test]
    fn chunk_count_respects_threads_and_minimum_size() {
        assert_eq!(chunk_count(10, 8), 1);
        assert_eq!(chunk_count(360_000, 3), 3);
        assert_eq!(chunk_count(360_000, 8), 3);
        assert_eq!(chunk_count(5_000_000, 4), 4);
    }

    #[test]
    #[should_panic]
    fn out_of_range_row_index_panics() {
        let mut out = vec![0.0; 2];
        row_sums_serial(&[1.0], &[5], &mut out, None);
    }
}
