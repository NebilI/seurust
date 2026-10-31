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
**Seurat 5.6.0**, then refreshed again after the Rust kernel optimizations in
`OPTIMIZATION_ROADMAP.md` (items 1–7): **34 of 38** kernels beat Seurat C++
on the published harness.

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
