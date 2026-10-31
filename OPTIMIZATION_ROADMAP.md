# seurust Rust optimization roadmap

Where the Rust kernels still lose to (or barely beat) Seurat 5.6.0's C++, why, and what to change. Sorted from most to least impactful.

**Status:** items **1–7 have landed** on `main`. The top-level [`README.md`](README.md) timing table is the post-opt remeasure (Seurat 5.6.0, `setThreads(3)`, 4-core machine): **34 of 38** kernels beat C++. Remaining losers: `CcaCrossprodMultiply` 0.71×, `EigenGramPCA` 0.78×, `SparseRowMeanVar` 0.80×, `fast_dist` 0.90×.

**Pre-opt baseline** (used for the projections below): the earlier 0.2.0 README table before items 1–7. Eight of 38 kernels were slower than C++ then: `LogNorm` 0.34×, `CcaCrossprodMultiply` 0.39×, `FastSparseRowScale` 0.41×, `FindWeightsC` 0.56×, `ComputeSNN` 0.59×, `Standardize` 0.85×, `RunModularityClusteringCpp` 0.94×, `RunModularityClusteringCpp_multi` 0.95×.

**How the estimates were made.** Items marked *measured* were prototyped in [`scripts/rust-prototype-bench/`](scripts/rust-prototype-bench/). That is a standalone Rust crate that runs the current kernel loop next to the proposed one, on the same input shapes as `scripts/bench-all-kernels.R` and on 4 cores. It does not link R, so a fresh untouched `malloc` stands in for R's `allocVector`, which also hands back unzeroed pages. That VM ran about 2× slower than the README machine overall, so only the **ratios** carry over. Items marked *estimated* come from reading the code. Projected speedups are relative to Seurat C++, computed as README Seurat time ÷ projected seurust time.

Run the prototype with:

```sh
cd scripts/rust-prototype-bench && cargo run --release
```

---

## Summary

| # | Optimization | Kernels | Current → projected vs C++ | Basis |
| --- | --- | --- | --- | --- |
| 1 | Write dense outputs straight into the R allocation | `FastSparseRowScale`, `Standardize`, `FastDenseRowScale`, `SCTPearsonResidualMatrix`, `EigenGramPCA`, `FindAllMarkers…` | `FastSparseRowScale` 0.41× → **~1.4×**; `Standardize` 0.85× → **~3×** | measured 3.7–4.8× |
| 2 | Column-parallel fill for `FastSparseRowScaleWithKnownStats` | `FastSparseRowScaleWithKnownStats` | 3.19× → **~14×** | measured 4.4–5.2× |
| 3 | Fused, parallel, copy-free `LogNorm` | `LogNorm` | 0.34× → **~1.6×** | measured 4.9–5.0× |
| 4 | `ComputeSNN`: prune before sorting, flat setup | `ComputeSNN` (+ `DirectSNNToFile`) | 0.59× → **~1.2×** | measured 1.3–1.6×; rest estimated |
| 5 | `FindWeightsC`: drop `String`/`HashMap`/`sprs`/`methods::new` | `FindWeightsC` | 0.56× → **~1.3×** | estimated 2–3× |
| 6 | Run `_multi` resolutions in parallel; trim the Louvain hot loop | `RunModularityClusteringCpp(_multi)` | `_multi` 0.95× → **~2×**; single 0.94× → ~1.05× | estimated |
| 7 | Serial fast path for tiny inputs | `CcaCrossprodMultiply` (+ other sub-0.1 ms kernels) | 0.39× → **~1×** | measured 2.4–2.9× |
| 8 | Stop copying `CscSlots` in and out; skip R-side `methods::new` | `RunUMISampling(PerCell)`, `GraphToNeighborHelper`, `ScoreHelper`, `SNN_SmallestNonzero_Dist` | ~1.1–1.6× → ~1.5–2.5× | estimated 1.2–1.6× |
| 9 | Cheaper row-statistics passes | `FastLogVMR`, `SparseRowVar`, `SparseRowMeanVar(Legacy)`, `FastExpMean`, `FindAllMarkers…` | e.g. `SparseRowMeanVar` 1.07× → ~1.4× | estimated 1.2–1.5× |
| 10 | Partial eigensolver / `dsyrk` Gram for `EigenGramPCA` | `EigenGramPCA` | 1.06× → ~1.5–2× | estimated |
| 11 | Parallel column copy in `FastRBind` | `FastRBind` | 1.28× → ~2× | estimated |
| — | Benchmark fairness fixes (not speedups) | `LogNorm`, all threaded kernels | changes the reported numbers | — |

Items 1–3 should come first. They are small, local edits to code that runs on every dataset (normalize, then scale). Measured, they recover 15–30 ms per call at the benchmark size, and that saving grows linearly with genes × cells. Items 1, 3 and 4 alone take three of the eight losing kernels past parity, and item 7 brings a fourth to about parity.

---

## 1. Write dense outputs straight into the R allocation (cross-cutting)

**Problem.** Many kernels build the result in a Rust `Vec`, then call `rmatrix_from_column_major` (`sparse.rs`), which allocates an R matrix and `memcpy`s the whole buffer into it on one thread:

```rust
let mut data = vec![0.0; n_sel * n_cells];           // 1st allocation, page-faulted while filling
data.par_chunks_mut(n_sel).for_each(/* fill */);
rmatrix_from_column_major(&data, n_sel, n_cells)       // 2nd allocation + serial 80 MB memcpy
```

For `FastSparseRowScale` (5,000 × 2,000) that is 80 MB faulted twice plus an 80 MB copy on one thread. The C++ version allocates once and fills in place. The fill itself is not the bottleneck.

`RMatrix::<f64>::new` already returns uninitialized R memory (`Robj::alloc_matrix`), so the fix is to allocate it first and run `par_chunks_mut` on `out.data_mut()`. Rayon workers only touch the raw slice and never call the R API, so this is safe.

A related double pass: `Doubles::new(n)` zero-fills element by element through an iterator, and `vec_to_doubles` / `into_r_dgcmatrix` / `log_norm` then `copy_from_slice` over it. Add one `alloc_doubles_uninit(n)` helper, using `Robj::alloc_vector(REALSXP, n)` the way `RMatrix::new` does, and use it everywhere a buffer is about to be fully overwritten.

**Where:** `fast_sparse_row_scale_impl`, `standardize_impl`, `fast_dense_row_scale_impl`, `sct_pearson_residual_matrix_impl`, `eigen_gram_pca_impl` (loadings/embeddings), `find_all_markers…` (`dimmed_matrix`), `graph_to_neighbor_helper_impl`, `vec_to_doubles`, `log_norm`.

**Measured** (prototype, current → direct write):

| Kernel | Current | Direct | Ratio |
| --- | ---: | ---: | ---: |
| `FastSparseRowScale` 5,000 × 2,000 | 37.3–42.3 ms | 9.6–11.4 ms | **3.7–3.9×** |
| `Standardize` 2,000 × 120 | 0.48 ms | 0.10–0.13 ms | **3.6–4.8×** |

**Projected:** `FastSparseRowScale` 25.7 ms → ~7 ms (0.41× → ~1.4×). `Standardize` 0.20 ms → ~0.05 ms (0.85× → ~3×). `FastDenseRowScale` and `SCTPearsonResidualMatrix` should improve by a similar factor (~2×).

**Extra for `FastSparseRowScale` / `FastDenseRowScale`:** the row-statistics pass that runs first is serial (and strided for the dense version). Parallelize it by **row range**: each thread owns a block of rows and scans all columns. That keeps the per-row summation order identical, so the current "exact" parity result survives. Splitting by columns and then reducing would change rounding. Worth roughly another 1 ms at benchmark size.

**Risk:** low. The values are unchanged; only where they are written changes.

## 2. Column-parallel fill for `FastSparseRowScaleWithKnownStats`

**Problem.** `fast_sparse_row_scale_with_known_stats_impl` parallelizes over **genes**. For each gene it writes all `n_cells` entries with a stride of `n_genes` (40 KB at benchmark size) into the column-major output. Every write touches a different cache line, and adjacent genes running on different threads share cache lines (false sharing). It also builds a full serial `RowIndex` first, then does the double allocation and copy from item 1.

**Fix:** use the same loop as `FastSparseRowScale`. Go parallel over **columns**: `copy_from_slice` a precomputed per-row "zero value" column, scatter that column's nonzeros, and write directly into the R matrix. `RowIndex` is no longer needed.

**Measured:** 46.3–56.4 ms → 8.9–12.7 ms, **4.4–5.2×**.
**Projected:** 28.5 ms → ~6 ms (3.19× → ~14× vs C++).
**Risk:** low. Each output element gets the same arithmetic, so results stay bit-identical.

## 3. Fused, parallel, copy-free `LogNorm`

**Problem.** Seurat's 5.6 kernel is threaded; ours is fully serial and makes extra passes:

1. `Doubles::new` zero-fills a new vector.
2. `copy_from_slice` copies `x` into it.
3. `col_sums()` scans all of `x` serially.
4. A second serial pass applies `log1p(x / colsum * scale)`.

**Fix:** allocate the output uninitialized and loop over columns in parallel. Each task sums its column, which is cache-hot, and immediately writes `log1p` into its own disjoint `p[c]..p[c+1]` range of the output. Each column is still summed sequentially, so results stay bit-identical.

**Measured:** 6.83–6.86 ms → 1.37–1.39 ms, **4.9–5.0×**.
**Projected:** 3.61 ms → ~0.75 ms (0.34× → ~1.6×).

Also worth doing: have the R wrapper skip `slot<-` (which runs `checkSlotAssignment` and may shallow-duplicate the object) by returning the finished `dgCMatrix` from Rust via `dgcmatrix_from_buffers`, reusing the input's `i`/`p` vectors.

**Risk:** low.

## 4. `ComputeSNN`: prune before sorting, flat setup

**Problem** (`compute_snn_accumulating_csc_parallel_from_data`):

- Every touched row of every column is **sorted before pruning**. With k = 20 that is about 400 rows per column. With `prune = 1/15` and realistic graphs, most of them are dropped right after the sort.
- `neighbors_per_cell_with_counts` builds `n_cells` separate heap `Vec`s with a linear duplicate search, all serial. In the prototype this setup alone took **0.51 ms of the 1.6 ms total**.
- Each column then allocates two more `Vec`s that are concatenated serially, and the result is copied again by `into_r_dgcmatrix` (`Doubles::new` zero-fill plus copy).
- `inv_cells` stores `usize` where `u32` would halve the memory traffic.

**Fix:**

1. `touched.retain(|r| scaled(acc[r]) >= prune)` **before** `sort_unstable`. The output is identical.
2. Replace the per-cell `Vec<Vec<(usize,u32)>>` with flat `n_cells × k` arrays. Dedupe each row's k neighbors with a tiny insertion sort or a stamp array.
3. Use `u32` indices.
4. Count nnz per column in a first parallel pass, then write `i`/`x` straight into R vectors (item 1 helper) instead of collecting per-column `Vec`s.

**Measured:** step 1 alone gives 2.49–2.52 ms → 1.57–1.91 ms (**1.3–1.6×**). Steps 2–4 should remove most of the remaining serial 0.5 ms (*estimated*).
**Projected:** 1.75 ms → ~0.8–0.9 ms (0.59× → ~1.2×). `DirectSNNToFile` gets the same speedup in its SNN half.
**Risk:** low. Steps 1–3 don't change the output. Step 4 only changes how the buffers are assembled.

## 5. `FindWeightsC`: drop `String`, `HashMap`, `sprs` and `methods::new`

**Problem** (`find_weights_impl` plus the R wrapper). For 400 cells × 600 anchors the actual arithmetic is tiny, so overhead dominates:

- `strings_to_str_vec` allocates a Rust `String` for every anchor name and every row name.
- Names are then grouped in a SipHash `HashMap<&str, Vec<usize>>`, and `cell_map` **clones** a `Vec` per anchor.
- The `min_dist == 0` path inserts every triplet into a second SipHash `HashMap<(usize,usize),f64>` for "last wins". It then goes through `sprs::TriMat` → `CsMat` (a `usize` conversion and a sort) → `CscSlots` → `to_r_list` (which **clones** `x`, `i` and `p`). Finally the R-side `CscFromList` calls `methods::new("dgCMatrix")`, which runs the full validity check.

**Fix:**

- Key names by CHARSXP pointer (`char_key`, already used in `row_merge_matrices_impl`) in a flat `Vec`, with no `String`s and no SipHash.
- Every triplet's column is its `cell`, so bucket triplets by cell and build CSC directly. For each column, sort that column's rows with a **stable** sort, keep the last duplicate, and normalize by the column sum in the same pass.
- Return the `dgCMatrix` from Rust with `build_dgcmatrix`.
- Optionally parallelize over `cells2` (the `nthreads` argument is currently ignored).

**Estimated:** 2–3×, 0.49 ms → ~0.2 ms (0.56× → ~1.3×).
**Risk:** low–medium. Last-wins duplicate semantics must match Eigen's `setFromTriplets` exactly. `FindWeightsC` is currently "exact" in the error table, so keep the parity test.

## 6. Modularity clustering: parallel resolutions, leaner hot loop

**`RunModularityClusteringCpp_multi`.** `run_modularity_clustering_multi_impl` loops over resolutions serially. Each resolution starts its own fresh `JavaRandom(random_seed)` and shares the network read-only through `Arc`, so the runs are fully independent. Running them with `par_iter` over `resolutions` gives **bit-identical** labels.

Estimated **~2–2.5×** for the benchmark's 3 resolutions on 3–4 cores: 0.95× → ~2×. In real use (`FindClusters` with a vector of resolutions) the gain scales with the number of resolutions.

**Single resolution** (`vos.rs`, `network.rs`). This is a line-for-line port, so the remaining gap is codegen:

- `run_local_moving_algorithm` has bounds checks on every `cluster[neighbor[k]]` / `edge_weight_per_cluster[l]` access. Use iterators over the `start..end` neighbor slice, or `get_unchecked` behind debug asserts.
- `calc_quality_function` re-sums every edge weight (`get_total_edge_weight`) on each iteration. Cache it in `Network`.
- `from_adjacency` clones `first_neighbor_index`, `neighbor` and `edge_weight` just to compute node weights. Compute them from the borrowed slices.
- The `i32` ↔ `usize` casts everywhere get in the way of vectorization. Store indices as `u32`/`usize` internally.

Estimated **5–15%**: 0.94× → ~1.0–1.05×. Low priority on its own, but cheap.

**Risk:** the `_multi` change is very low risk. The hot-loop edits must not change iteration order, or RNG parity with C++ is lost.

## 7. Serial fast path for tiny inputs

**Problem.** `CcaCrossprodMultiply` (80 × 200) takes 0.01 ms in C++. Ours spends about as long dispatching onto the rayon pool as computing. Other kernels that always go through rayon have the same issue at small sizes: `CountAnchorSharedNeighbors`, `align_loading_signs` (which runs `par_iter_mut` over 120 values), and `fast_cov_mats` scaling.

**Fix:** a shared `if work < THRESHOLD { serial } else { rayon }` helper, or `with_min_len` on the parallel iterators. A threshold of roughly 50–100k flops is a sensible starting point.

**Measured** (`CcaCrossprodMultiply`): 0.020–0.021 ms → 0.007–0.009 ms, **2.4–2.9×**. **Projected:** 0.39× → ~1×. The absolute saving is microseconds, but a "slower than C++" row in the README is visible.

## 8. Stop copying `CscSlots` in and out; skip R-side `methods::new`

**Problem.** `CscSlots::from_r` copies `x`, `i` and `p` element by element through iterators. Several kernels then copy again:

- `RunUMISampling` / `RunUMISamplingPerCell`: `from_r` copies → mutate → `to_r_list` (clones all three vectors, then `Doubles::new` zero-fill plus copy) → R `CscFromList` → `methods::new` validity scan. That is about 5 full passes over the matrix in addition to the actual sampling.
- `GraphToNeighborHelper`: `from_r` → `to_cs_mat` (another copy, widening to `usize`) → `to_csr` → a per-row `Vec<usize>` order allocation.
- `ScoreHelper`: `from_r` → `to_cs_mat` → per-column `Vec`s.
- `SNN_SmallestNonzero_Dist`: `from_r`, plus `ndarray_from_rmatrix` copies the embedding matrix.

**Fix:**

- Read through a `CscView` (zero-copy) and write the output into fresh R vectors. For UMI sampling, only `x` changes: allocate `x_out` and reuse the input `i`/`p` R vectors directly in `dgcmatrix_from_buffers`.
- Do the CSC→CSR transpose with a counting-sort pass (as `RowIndex` already does) instead of going through `sprs`.
- Reuse one scratch buffer per rayon task instead of per-row `Vec`s.

**Estimated:** 1.2–1.6× each. UMI sampling is bounded by the inherent per-entry `Rf_runif` call, which must stay serial for RNG parity.

## 9. Cheaper row-statistics passes

**Problem** (`par_row_sums`, `par_row_sums_and_counts`, `sparse_row_mean_var_impl` and friends):

- `fold` creates a fresh `n_genes`-sized accumulator, or two of them, for **every rayon split**, not every thread. Each one must be zeroed and then reduced. At 5,000 genes × 2,000 columns that is a noticeable share of a 0.25 ms kernel.
- `FastLogVMR` and `SparseRowVar` scan the matrix twice, and `FastLogVMR` recomputes `expm1` on every nonzero in the second pass.

**Fix:**

- Split columns into exactly `current_num_threads()` contiguous chunks. That gives one accumulator per thread, as the hand-rolled split in `integrate_data_impl` already does. A deterministic chunking also makes the reduction order reproducible across runs.
- In `FastLogVMR`, store the `expm1` values from pass 1 in a reusable buffer, or interleave two-pass per column chunk, so pass 2 reads them back.
- Fuse the count and sum passes where the C++ formula allows it without changing rounding.

**Estimated:** 1.2–1.5× (e.g. `SparseRowMeanVar` 1.07× → ~1.4×, `FastLogVMR` 4.8× → ~6×).
**Risk:** low, but the reduction order affects the last bits of the results, so rerun `SEURUST_BENCH_MODE=error`.

## 10. `EigenGramPCA`: partial eigensolver and symmetric Gram

**Problem.** `nalgebra::SymmetricEigen` computes **all** 120 eigenpairs, and the full result is then sorted. Seurat uses Spectra to compute only `npcs`. The Gram matrix uses a general `dgemm` (twice the flops of a symmetric rank-k update), and the outputs are copied (item 1).

**Fix:** call R's LAPACK `dsyevr` with `RANGE='I'` (top-`npcs` eigenpairs only). It is available to any R package via `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` in `Makevars`. Alternatively, implement a small Lanczos iteration with restarts. Compute the Gram with `dsyrk`.

**Estimated:** 1.5–2× (1.06× → ~1.5–2×). The gap widens on realistic sizes: 2,000 variable features and 50 PCs is where a full eigendecomposition really hurts.
**Risk:** medium. It changes numerics. `dsyevr` should *improve* agreement with Spectra (currently 7.2 significant figures), but verify it.

## 11. Parallel column copy in `FastRBind`

`fast_rbind_impl` copies 9.6 MB on a single thread. Use `par_chunks_mut(n)` over the output columns, with one task per column block. That is memory-bound, so expect about **1.5×** on 3–4 cores (1.28× → ~2×). Trivial and low risk.

---

## Benchmark fairness fixes (not speedups)

These change what the README reports without changing the kernels:

- **`LogNorm` is not compared like for like.** The harness calls Seurat's raw `LogNorm(sparse@x, sparse@p, …)`, which returns a numeric vector. On the seurust side it calls the full R wrapper (`seurust::LogNorm(sparse, …)@x`), which also does the class check, `slot<-` and S4 handling. Either time `seurust:::log_norm` on the slots, or wrap the Seurat call the same way.
- **Thread counts differ.** Seurat runs with `setThreads(3)`, while seurust ignores every `nthreads` argument (`_nthreads`) and uses rayon's global pool, which defaults to all 4 cores. Honor `nthreads` by running inside a cached `rayon::ThreadPool` of that size (e.g. `pool.install(|| …)`, keyed by thread count). That makes the comparison fair, and it also respects users who limit threads on shared HPC nodes.
- **Sub-0.1 ms kernels use 5–7 repeats.** That is too few to resolve a 0.01 ms difference. Use `times = 100` or more for kernels under 0.1 ms so rows like `CcaCrossprodMultiply` and `CountAnchorSharedNeighbors` are stable from run to run.
- **One input size per kernel.** Add a "large" tier, for example 20,000 × 50,000 at about 7% density, for the preprocessing kernels. Items 1–3 and 9 scale linearly with matrix size, so their real-world payoff is much bigger than the 5,000 × 2,000 numbers suggest.

## Correctness guardrails for all of the above

- Keep the bit-for-bit "exact" rows exact. Prefer parallel splits that keep each output element's accumulation order: row-range splits for row statistics, and per-column work for column outputs.
- After every change, rerun `SEURUST_BENCH_MODE=error Rscript scripts/bench-all-kernels.R` and the `testthat` parity suite (`seurust/tests/testthat/test-seurat56.R`).
- Writing to R memory from rayon workers is fine as long as workers touch only raw slices. Every R API call (allocation, attributes, `Rf_runif`) must stay on the main thread.
