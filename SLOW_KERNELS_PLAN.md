# Plan: the four kernels still slower than Seurat C++

After roadmap items 1–7 ([`OPTIMIZATION_ROADMAP.md`](OPTIMIZATION_ROADMAP.md)), 34 of 38 kernels beat Seurat 5.6.0. This plan covers the four that still lose, in the order we intend to fix them. Each step is a separate change that can be reviewed, benchmarked and merged on its own.

| Step | Kernel | Today vs C++ | Projected vs C++ | Basis | Risk |
| --- | --- | ---: | ---: | --- | --- |
| 1 | `SparseRowMeanVar` | 0.80× | **~2.8×** (2.5–3.0×) | measured 3.5–3.9× | low, bit-identical on the benchmark |
| 2 | `EigenGramPCA` | 0.78× | **~2.7×** (2.4–3.0×) | Gram measured 3.1×; eigensolver estimated from a measured lower bound | medium, changes numerics slightly |
| 3 | `fast_dist` | 0.90× | **~1.3×** (1.2–1.4×) | arithmetic measured; R object overhead modelled | low, bit-identical |
| 4 | `CcaCrossprodMultiply` | 0.71× | **~0.95×** (0.9–1.05×) | compute measured 1.25–1.4×; R wrapper cost estimated | low, bit-identical |
| 5 (optional) | Roll step 1's accumulator out to the other row-statistics kernels | already 1.3–7× | another ~1.5–3× each | same mechanism as step 1 | low |

If steps 1–3 land as projected, 37 of 38 kernels beat C++. `CcaCrossprodMultiply` is a 10 µs call, so step 4 is about reaching parity rather than saving any real time.

## How the numbers were made

The prototypes are in [`scripts/rust-prototype-bench/src/round2.rs`](scripts/rust-prototype-bench/src/round2.rs). They use the same input shapes as `scripts/bench-all-kernels.R` and run on 4 cores:

```sh
cd scripts/rust-prototype-bench && cargo run --release -- round2
```

The prototype does not link R. That VM had no R or Seurat installed and runs about 1.9× slower than the README machine. For example, the current `EigenGramPCA` pipeline took 1.60 ms there and 0.81 ms in the README. So only the **ratios** between current and proposed code carry over. Projections are the README seurust time divided by the measured ratio, then compared to the README Seurat time. Each step must be re-benchmarked against Seurat with `scripts/bench-all-kernels.R` before we claim the new number.

---

## Step 1: `SparseRowMeanVar` — one accumulator per thread

**Why it is slow.** `sparse_row_mean_var_impl` uses rayon `fold` / `reduce`. `fold` makes a fresh pair of `n_genes`-long accumulators (2 × 40 KB at 5,000 genes) for **every rayon split**, not every thread. Each one is zeroed, then all of them are merged pairwise. At this size that bookkeeping costs more than the actual scatter-add: the current parallel code is **slower than a plain serial loop**.

| Variant (prototype, 5,000 × 2,000, 5% nonzero) | Time |
| --- | ---: |
| Current `fold` / `reduce` | 0.46–0.51 ms |
| Plain serial loop | 0.32 ms |
| **Fixed chunks, 3-thread pool** (Seurat runs with `setThreads(3)`) | **0.12–0.14 ms** |
| Fixed chunks, 4 threads | 0.11–0.17 ms |

**Change.**

1. Split the columns into exactly `current_num_threads()` contiguous chunks, balanced by nonzero count using `p`.
2. Give each chunk one interleaved `[sum, sumsq]` accumulator, so each nonzero updates one cache line instead of two.
3. Reduce in parallel by **gene range**: each thread owns a slice of genes and adds up all chunk accumulators for those genes, then writes the mean and variance directly.
4. Write `mean` and `variance` straight into R vectors (`doubles_uninit`) instead of building `Vec`s and copying them with `Doubles::from_values`.

**Projected.** 0.26 ms → ~0.07–0.08 ms, so 0.80× → **~2.8×** vs C++.

**Correctness.** On the benchmark input, the counts are integers, so the sums are exact in any order and the output is bit-identical (the prototype checks this). On non-integer data the last bits can change, but the chunking is deterministic for a given thread count. The current `fold` split points depend on work stealing, so its rounding is not reproducible across runs. The error table should keep `SparseRowMeanVar` at "exact".

**Files.** `seurust/src/rust/src/data_manipulation/sparse_fns.rs` (new shared helper plus `sparse_row_mean_var_impl`), `seurust/src/rust/src/lib.rs` (`sparse_row_mean_var` output).

## Step 2: `EigenGramPCA` — partial eigensolver and parallel Gram

**Why it is slow.** In the prototype, 120 features × 400 cells × 10 PCs:

| Phase | Current | Share |
| --- | ---: | ---: |
| Gram `X Xᵀ` (serial `matrixmultiply::dgemm`) | 0.235 ms | 15% |
| `nalgebra::SymmetricEigen`: **all 120** eigenpairs, then sort | 1.32 ms | **83%** |
| Embeddings `Xᵀ U` | 0.032 ms | 2% |

Seurat only computes the top `npcs` eigenpairs (Spectra). We compute all of them with nalgebra's unblocked QR iteration, which also accumulates every eigenvector.

### 2a. Parallel, symmetric Gram

Split the Gram into column blocks on the rayon pool and compute only the lower block triangle, then mirror it. This also starts honoring `nthreads`.

| Gram variant | 120 × 400 | 1,000 × 3,000 |
| --- | ---: | ---: |
| Current serial `dgemm` | 0.235 ms | 123–136 ms |
| 4 column blocks | 0.083 ms | 32–46 ms |
| **Lower block triangle, parallel** | **0.075 ms** | **24–26 ms** |

The result is bit-identical to the current Gram in both cases, because each output element keeps the same inner-product order inside `matrixmultiply`.

### 2b. Top-`npcs` eigenpairs only

Call R's own LAPACK `dsyevr` with `RANGE = 'I'` for eigenpairs `n − npcs + 1 … n`. R already uses `dsyevr` for `eigen(symmetric = TRUE)`, so it is always available. Link it by adding `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` to `PKG_LIBS` in `Makevars.in` / `Makevars.win.in` and declaring the Fortran symbol in Rust. `dsyevr` reduces to tridiagonal form once (O(n³), blocked), then finds only the requested eigenpairs (O(n·k)) and back-transforms them (O(n²·k)).

The prototype cannot call LAPACK without R. Instead it timed nalgebra's Householder tridiagonalization alone, the same O(n³) step `dsyevr` has to do: **0.28 ms vs 1.32 ms** for the full solve. Adding the small subset solve and back-transform gives an estimated **0.30–0.40 ms** for the eigen phase.

**Not recommended: a hand-written Lanczos solver.** On the benchmark's `rnorm` data the top eigenvalues are tightly clustered. Plain Lanczos needed 60 steps to reach 1e-5 residuals at 120 × 120, which took 0.68 ms. At 1,000 × 1,000 it did not converge at all without restarts. Matching Spectra's restarted Lanczos is a lot of code for less benefit than `dsyevr`.

**Projected.** In the prototype, the pipeline goes from 1.60 ms to about 0.08 + 0.35 + 0.03 ≈ 0.46 ms, which is about 3.5× faster. In the README, 0.81 ms → ~0.23 ms (0.21–0.26 ms), so 0.78× → **~2.7×** vs C++. 2a alone only gets to ~0.87×, so 2b is the step that matters.

**At realistic sizes the gain is bigger.** With 1,000 features × 3,000 cells × 30 PCs, the current pipeline takes 0.86–0.90 s, of which 0.72–0.77 s is the full eigensolve. With 2a and 2b it should take ~0.25–0.3 s, about 3× faster. Seurat's default of 2,000 variable features would widen the gap further, because the full solve is O(n³) in features.

**Correctness.** `dsyevr` is accurate to machine precision, so our side becomes more accurate. The current 7.2-significant-figure agreement is limited by Spectra's convergence tolerance, so the error-table number may not move much. Keep the sign alignment (`align_loading_signs`) and the descending order. Check the result against `eigen(crossprod(t(X)), symmetric = TRUE)` in a `testthat` test.

**Files.** `seurust/src/rust/src/pca.rs`, `seurust/src/Makevars.in`, `seurust/src/Makevars.win.in`.

## Step 3: `fast_dist` — build the result list without extendr bookkeeping

**Why it is slow.** The benchmark returns a list of 800 numeric vectors of length 10. The arithmetic is small: 8,000 distances × 20 dimensions. Most of the time goes to creating R objects:

- extendr 0.9 tracks every `Robj` it creates or reads in a global `Mutex<HashMap<SEXP, …>>` (`ownership::protect` / `unprotect`, SipHash). `fast_dist` goes through it about 1,600 times per call: once for each list element read through `n.values()`, once for each `Doubles::new`, and again while building the result `List`.
- When names exist, they are converted to Rust `&str`s and then back into CHARSXPs by `set_names`.

| Prototype measurement (800 × 20, 10 neighbors) | Time |
| --- | ---: |
| Distance arithmetic, current strided loop | 0.074 ms |
| Same, transposed to row-major first | 0.083 ms (slower: both inputs already fit in L2 cache) |
| Same, 4 neighbors per pass | 0.069 ms (6% faster, bit-identical) |
| Modelled extendr protect + unprotect for ~1,600 objects | **0.070 ms** |

**Change.**

1. Allocate the output `VECSXP` once, protect it once, and for each cell `Rf_allocVector(REALSXP, len)` then `SET_VECTOR_ELT` straight into it. Read neighbor vectors with `VECTOR_ELT` and `INTEGER` / `REAL` on raw `SEXP`s. With this, no per-element `Robj` exists, so extendr's ownership table is never touched.
2. Copy the input's `names` attribute pointer onto the output instead of round-tripping through `&str`.
3. Compute 4 neighbors per pass over the dimensions. It is cheap and bit-identical.
4. Fold the `ngraph_size >= 4096` branch into the same code. Compute distances into one flat buffer in parallel when the work is large, then do the serial R allocation and `memcpy` per cell.

**Projected.** In the README, 0.19 ms → ~0.13–0.15 ms, so 0.90× → **~1.3×** vs C++. This is the least certain projection because the extendr overhead is modelled rather than measured inside R, so measure it in R first.

**Bonus.** The same raw-`SEXP` list building applies wherever we return many small R objects, for example the per-cell lists in `GraphToNeighborHelper`.

**Files.** `seurust/src/rust/src/fast_nn_dist.rs`, possibly `seurust/R/native.R` (`as.list` is a no-op for lists, so it can stay).

## Step 4: `CcaCrossprodMultiply` — shorter dependency chains and a leaner wrapper

**Why it is slow.** Both sides take about 10 µs. Ours is roughly half compute and half call overhead:

- The second loop computes 200 dot products of length 80. Each is a single serial chain of dependent adds, so the CPU stalls on add latency.
- The R wrapper calls `AsDoubleMatrix` twice and `as.double` once. Each is an R closure call doing `is.matrix` / `is.double` checks. That costs a few µs, which is a large fraction of a 10 µs call.

| Prototype measurement | Time |
| --- | ---: |
| Current compute | 6.8–7.7 µs |
| **4 output columns per pass** (4 independent chains, each still summed in row order) | **5.5 µs**, bit-identical |

**Change.**

1. Compute four output columns per pass, as prototyped.
2. In `CcaCrossprodMultiply`, skip the coercion helpers when the inputs are already double matrices (`if (is.double(left) && is.matrix(left))`, etc.), or do the type check in Rust and only fall back to R coercion when it fails.
3. Benchmark fairness: raise this case from 7 to ≥ 200 repeats in `bench-all-kernels.R`. At 7 repeats a 1–2 µs difference is noise.

**Projected.** About 10 µs → 7–8 µs, so 0.71× → **~0.95×** (0.9–1.05×). This is near parity and could land either side of 1.0×. It saves only a few µs per call, so it goes last.

**Files.** `seurust/src/rust/src/integration.rs`, `seurust/R/native.R`, `scripts/bench-all-kernels.R`.

## Step 5 (optional): reuse step 1's accumulator everywhere

The same per-split `fold` pattern lives in `par_row_sums`, `par_row_sums_and_counts`, `sparse_row_var2_impl`, `sparse_row_var_std_legacy_rows_impl`, `find_all_markers_sparse_fold_change_stats_impl` and SCT `residual_stats`. These already beat C++, but they back `FastExpMean`, `SparseRowVar`, `FastLogVMR`, `SparseRowVarStd`, `SparseRowMeanVarLegacy`, `SparseRowVarStdLegacyRows`, `FindAllMarkersSparseFoldChangeStats` and `SCTResidualStatsAndCorrected`. Switching them to the step-1 helper should give a similar ~1.5–3× on the row-statistics part of each kernel. This is the same as roadmap item 9.

---

## Verification for every step

1. `cargo test` in `seurust/src/rust`, plus a new unit test that compares the new code path with the old one (bit-for-bit where the step promises it).
2. `SEURUST_BENCH_MODE=error Rscript scripts/bench-all-kernels.R`. "exact" rows must stay exact, and `EigenGramPCA` must not lose significant figures.
3. `Rscript scripts/bench-all-kernels.R` on the README machine for the before/after timing table.
4. `testthat` parity suite: `seurust/tests/testthat/test-seurat56.R`.

The VM these prototypes ran on has no R or Seurat, so steps 2 and 3 need an R + Seurat 5.6.0 environment to verify. The Docker images in [`docker/`](docker/) are the easiest way to get one.
