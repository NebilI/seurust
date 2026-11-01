# seurust (unreleased)

Further kernel work from [`SLOW_KERNELS_PLAN.md`](https://github.com/NebilI/seurust/blob/main/SLOW_KERNELS_PLAN.md) for the four functions that still lost to Seurat 5.6.0 C++ after 0.3.0.

## Improvements

* **`SparseRowMeanVar` and other row-statistics kernels** — Columns are split into one nnz-balanced chunk per rayon thread. Each thread owns a single accumulator (no per-split `fold` allocations). `SparseRowMeanVar` writes `mean` / `variance` straight into uninitialized R vectors. The same helper is used by `FastExpMean`, `SparseRowVar`, `SparseRowVar2`, `SparseRowVarStd`, `FastLogVMR`, `SparseRowMeanVarLegacy`, `SparseRowVarStdLegacyRows`, `FindAllMarkersSparseFoldChangeStats`, and `SCTResidualStatsAndCorrected`.
* **`EigenGramPCA`** — Parallel lower-triangle Gram (`X Xᵀ`) honoring `nthreads`, then R's LAPACK `dsyevr` for only the top `npcs` eigenpairs instead of a full nalgebra solve.
* **`fast_dist`** — Builds the result `VECSXP` with the R C API (no per-element extendr ownership table), copies the input `names` attribute pointer, and computes four neighbor distances per pass over the dimensions.
* **`CcaCrossprodMultiply`** — Four output columns per pass (bit-identical) and a wrapper that skips `AsDoubleMatrix` / `as.double` when the inputs are already double matrices or a double vector.

# seurust 0.3.0

Performance release: seven kernel optimizations from
[`OPTIMIZATION_ROADMAP.md`](https://github.com/NebilI/seurust/blob/main/OPTIMIZATION_ROADMAP.md),
plus a single publish path for GitHub Releases, crates.io, and r-universe.
Public R signatures and per-entry arithmetic are unchanged from **0.2.0**
(results stay bit-identical or within the same numeric-error bands as before).

## Improvements

* **In-place dense writes** — Large dense outputs
  (`Standardize`, `FastSparseRowScale`, `FastDenseRowScale`,
  `FastSparseRowScaleWithKnownStats`, `SCTPearsonResidualMatrix`,
  `EigenGramPCA`, `GraphToNeighborHelper`, `LogNorm`, and related helpers)
  fill uninitialized R allocations directly instead of building a Rust `Vec`
  and copying it on one thread.
* **Column-parallel `FastSparseRowScaleWithKnownStats`** — Precomputes per-row
  means, inverse sdevs, and zero-column values, then fills each output column
  in parallel from CSC nonzeros (drops the serial `RowIndex` and gene-strided
  writes).
* **Fused, parallel `LogNorm`** — Sums each CSC column and writes
  `log1p(x / colsum * scale)` in the same pass, with columns on rayon and an
  uninitialized output buffer (no extra copy or separate `col_sums` scan).
* **Faster `ComputeSNN`** — Builds a flat `n × k` neighbor index with `u32`
  ids, inverts it once, prunes SNN edges before the per-column sort, and writes
  integer CSC slots into uninitialized R allocations.
* **Leaner `FindWeightsC`** — Keys integration rownames by interned `CHARSXP`,
  buckets weights by cell, applies Eigen last-wins with a stable per-column
  sort, and returns a `dgCMatrix` from Rust (no `String` / `HashMap` / `sprs`
  path; R wrapper skips `methods::new`).
* **Parallel multi-resolution modularity** —
  `RunModularityClusteringCpp_multi` runs independent resolutions with rayon
  (each still starts its own `JavaRandom`, so labels stay bit-identical). The
  Louvain local-moving loop caches total edge weight, derives node weights from
  borrowed slices, and walks neighbor ranges with slice iterators.
* **Serial fast path for tiny kernels** — Shared helper skips rayon when work
  is under a 64k-op floor. At benchmark size that covers
  `CcaCrossprodMultiply`, `CountAnchorSharedNeighbors`, `FastCovMats` scaling,
  and loading-sign flips; large kernels stay parallel.
* **Optimization roadmap** — Added `OPTIMIZATION_ROADMAP.md` and
  `scripts/rust-prototype-bench/` so remaining kernel work is measured against
  the seurust 0.2.0 / Seurat 5.6.0 baseline.

## Fixes / packaging

* **One release workflow** — Dropped the CRAN submit path and separate
  crates.io / r-universe Actions. **Release seurust** now creates the GitHub
  Release, publishes the Rust crate, builds the R tarball, and upserts
  r-universe in one job graph.
* **Release notes from NEWS** — GitHub Release bodies are generated from the
  matching `# seurust <version>` section in this file (see
  [`RELEASING.md`](https://github.com/NebilI/seurust/blob/main/RELEASING.md)
  and `scripts/prepare-release.sh`).

## Compatibility

Requires **Seurat (>= 5.6.0)** for the same kernel surface as 0.2.0.
`AnnoySearchCpp` remains unported. `WriteEdgeFile`, `DirectSNNToFile`, and
`SparseRowVar2` remain available in seurust even though Seurat 5.6.0 dropped
those R exports.

# seurust 0.2.0

Ported the native kernels added in **Seurat 5.6.0** and refreshed parity tests
and kernel benchmarks against that release.

## New kernels (Seurat 5.6.0)

* `FastDenseRowScale` — dense row center/scale/clip used by `ScaleData`.
* `SparseRowMeanVar` / `SparseRowMeanVarLegacy` / `SparseRowVarStdLegacyRows`
  — shared sparse VST statistics.
* `FindAllMarkersSparseFoldChangeStats` — bulk fold-change / detection prep.
* `EigenGramPCA` — Gram-matrix PCA backend used by default `RunPCA`.
* `SCTResidualStatsAndCorrected` / `SCTPearsonResidualMatrix` — SCTransform
  Pearson residuals and corrected counts.
* `CcaCrossprodMultiply` / `CountAnchorSharedNeighbors` — CCA integration
  matvec and threaded anchor scoring.
* `RunModularityClusteringCpp_multi` — reuse one SNN network across
  resolutions.

## Signature updates

Existing wrappers accept the extra `nthreads` (and, where Seurat added it,
`features`) arguments from 5.6.0. Convenience matrix-taking APIs such as
`LogNorm()` and `FastSparseRowScale()` keep working with seurust 0.1.0
positional calls. Seurat 5.6.0's C++ `LogNorm` now returns only the `x` slot;
`seurust::LogNorm` still returns a full `dgCMatrix`, matching
`Seurat::LogNormalize`.

`AnnoySearchCpp` is **not** ported (on-disk Annoy index search stays in Seurat).
`WriteEdgeFile`, `DirectSNNToFile`, and `SparseRowVar2` remain available in
seurust even though Seurat 5.6.0 dropped those R exports.

Benchmark tables in the repository `README.md` were remeasured against
**Seurat 5.6.0**.

# seurust 0.1.0

Fresh baseline release after syncing this fork with upstream **Seurat 5.5.1** and
resetting the seurust version line to **0.1.0** (replacing prior 0.1.x GitHub
release tags).

## Highlights

* **30** native kernels ported to Rust with parity tests against Seurat's C++
  implementations.
* Benchmark tables in the repository `README.md` were refreshed against
  **Seurat 5.5.1**. Per-kernel medians can be reproduced with
  `scripts/bench-all-kernels.R` against an installed copy of Seurat.
* Workflow-level timings (PBMC tutorials and a standard analysis script) used
  the same Seurat 5.5.1 and seurust 0.1.0 pairing. Those helper scripts are
  not part of the installed package.

## Notes

* crates.io may still list older 0.1.x crate versions that were published before
  this reset; they cannot be overwritten. New publishes use this **0.1.0**
  crate version when the registry allows it.
