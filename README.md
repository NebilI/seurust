# seurust

**Rust native kernels for [Seurat](https://satijalab.org/seurat) — same function signatures, less time in the hot path.**

[seurust](seurust/) is an R package that reimplements Seurat's performance-critical native routines with Rust and [extendr](https://extendr.github.io/). Install it next to Seurat, keep your analysis code, and call the same kernel names.

> **Drop-in by design.** seurust exposes the same function signatures as Seurat's internal C++ layer (`LogNorm`, `FastSparseRowScale`, `ComputeSNN`, `IntegrateDataC`, `EigenGramPCA`, and more). Parity checks compare results with the installed Seurat package ([`seurust/tests/testthat/`](seurust/tests/testthat/), [`scripts/bench-all-kernels.R`](scripts/bench-all-kernels.R)).

[![seurust CI](https://github.com/NebilI/seurust/actions/workflows/seurust_checks.yaml/badge.svg)](https://github.com/NebilI/seurust/actions/workflows/seurust_checks.yaml)
[![r-universe](https://NebilI.r-universe.dev/badges/seurust)](https://NebilI.r-universe.dev/seurust)

Publishing notes: [`r-universe/README.md`](r-universe/README.md).
PR tests: `seurust Checks`. Releases: Actions → **Release seurust**.

---

## Why use seurust?

Single-cell pipelines spend a surprising amount of time in a handful of native kernels: log-normalization, variable-feature statistics, scaling sparse matrices, building shared-nearest-neighbor graphs, and batch-integration weighting. Those routines dominate preprocessing and graph construction on large datasets.

seurust targets exactly that layer:

- **Same Seurat kernel names** — no new object model
- **Validated parity** — tests against the installed Seurat package on ported kernels
- **Measurable speedups** — Rust wins on the preprocessing kernels that run on every dataset (see benchmarks below)
- **Open development** — install from GitHub or r-universe

Seurat remains the user-facing analysis package. **seurust** is the native backend you install alongside it.

---

## Layout

| What | Where |
|------|--------|
| R package | [`seurust/`](seurust/) |
| Rust crate | [`seurust/src/rust/`](seurust/src/rust/) |
| Offline Rust crates | [`seurust/src/rust/vendor.tar.xz`](seurust/src/rust/vendor.tar.xz) |
| Kernel benchmark | [`scripts/bench-all-kernels.R`](scripts/bench-all-kernels.R) |
| Release process | [`RELEASING.md`](RELEASING.md) |
| Release helpers | [`scripts/prepare-release.sh`](scripts/prepare-release.sh), [`scripts/extract-release-notes.sh`](scripts/extract-release-notes.sh) |
| Release-tag cleanup | [`scripts/prune-seurust-releases.sh`](scripts/prune-seurust-releases.sh) |
| Docker development images | [`docker/`](docker/) |
| r-universe registry notes | [`r-universe/`](r-universe/) |

Seurat itself is installed separately (this repository does not contain Seurat's sources).

---

## Performance vs Seurat (C++)

Median runtime for every ported native function, Seurat's C++ kernels against seurust. **Speedup = Seurat time ÷ seurust time.** Values **above 1.0× mean seurust is faster**.

### seurust 0.3.0 + unreleased #45 kernels (`f1443e1`, Seurat 5.6.0 kernels)

Measured against **Seurat 5.6.0** with `setThreads(3)` on the machine described in [Machine](#machine), using seurust at commit `f1443e1` (`main` after [#45](https://github.com/NebilI/seurust/pull/45), which is not yet in a release). Seurat 5.6.0 dropped `SparseRowVar2`, `WriteEdgeFile`, and `DirectSNNToFile` from its R exports (seurust still implements them, but they are not in this comparison). `AnnoySearchCpp` is not ported. **36 of 38** kernels were faster. The four that trailed Seurat before #45 are now ahead: `SparseRowMeanVar` (1.55×), `CcaCrossprodMultiply` (1.52×), `fast_dist` (1.29×), and `EigenGramPCA` (1.28×). The two still slower are sub-millisecond: `row_sum_dgcmatrix` (0.92×) and `row_mean_dgcmatrix` (0.98×). Kernels within about ±10% of 1.00× can change sides between runs (for example `RunModularityClusteringCpp` at 1.07×), so treat those as parity.

| Function | Problem size | Seurat C++ | seurust Rust | Speedup |
| --- | --- | ---: | ---: | ---: |
| `LogNorm` | 5,000 × 2,000 sparse, 5% nonzero | 1.24 ms | 0.86 ms | **1.43×** |
| `Standardize` | 2,000 × 120 dense | 0.27 ms | 0.10 ms | **2.74×** |
| `FastCov` | 2,000 × 80 dense | 0.93 ms | 0.36 ms | **2.57×** |
| `FastCovMats` | 2,000 × 80 and 2,000 × 40 dense | 0.75 ms | 0.45 ms | **1.67×** |
| `FastRBind` | two 4,000 × 150 dense | 2.72 ms | 2.20 ms | **1.24×** |
| `RowVar` | 2,000 × 120 dense | 0.13 ms | 0.07 ms | **1.82×** |
| `FastExpMean` | 5,000 × 2,000 sparse, 5% nonzero | 4.49 ms | 0.76 ms | **5.91×** |
| `SparseRowVar` | 5,000 × 2,000 sparse, 5% nonzero | 2.66 ms | 0.50 ms | **5.30×** |
| `SparseRowVarStd` | 5,000 × 2,000 sparse, 5% nonzero | 1.00 ms | 0.41 ms | **2.41×** |
| `SparseRowMeanVar` | 5,000 × 2,000 sparse, 5% nonzero | 0.42 ms | 0.27 ms | **1.55×** |
| `SparseRowMeanVarLegacy` | 5,000 × 2,000 sparse, 5% nonzero | 0.57 ms | 0.40 ms | **1.40×** |
| `SparseRowVarStdLegacyRows` | 5,000 × 2,000 sparse, 80 rows | 0.56 ms | 0.23 ms | **2.45×** |
| `FastLogVMR` | 5,000 × 2,000 sparse, 5% nonzero | 7.42 ms | 2.11 ms | **3.52×** |
| `FastSparseRowScale` | 5,000 × 2,000 sparse, 5% nonzero | 14.0 ms | 7.08 ms | **1.97×** |
| `FastDenseRowScale` | 800 × 300 dense | 0.49 ms | 0.28 ms | **1.74×** |
| `FastSparseRowScaleWithKnownStats` | 5,000 × 2,000 sparse, 5% nonzero | 104 ms | 5.73 ms | **18.20×** |
| `RowMergeMatrices` | 1,500 × 400 and 1,200 × 350 CSR | 0.77 ms | 0.27 ms | **2.86×** |
| `ReplaceColsC` | 2,000 × 600 sparse, replace 20 columns | 0.30 ms | 0.12 ms | **2.59×** |
| `GraphToNeighborHelper` | 2,000 cells, 30 neighbors, symmetric | 0.79 ms | 0.66 ms | **1.21×** |
| `RunUMISampling` | 2,000 × 1,500 counts, 8% nonzero | 3.87 ms | 3.02 ms | **1.28×** |
| `RunUMISamplingPerCell` | 2,000 × 1,500 counts, 8% nonzero | 3.25 ms | 2.94 ms | **1.11×** |
| `ComputeSNN` | 1,500 cells, k = 20 | 1.30 ms | 1.04 ms | **1.25×** |
| `IntegrateDataC` | 800 genes × 200 cells, 60 anchors | 0.86 ms | 0.59 ms | **1.45×** |
| `FindWeightsC` | 400 cells, 600 anchors, k = 10 | 0.40 ms | 0.16 ms | **2.42×** |
| `ScoreHelper` | 300 cells × 16 dimensions | 0.63 ms | 0.51 ms | **1.22×** |
| `SNN_SmallestNonzero_Dist` | 400 cells × 10 dimensions | 0.41 ms | 0.33 ms | **1.23×** |
| `RunModularityClusteringCpp` | 250-cell SNN, 10 iterations | 4.93 ms | 4.60 ms | **1.07×** |
| `RunModularityClusteringCpp_multi` | 250-cell SNN, 3 resolutions | 2.57 ms | 2.55 ms | **1.01×** |
| `FindAllMarkersSparseFoldChangeStats` | 5,000 × 2,000 sparse, 4 groups | 1.93 ms | 0.82 ms | **2.36×** |
| `CcaCrossprodMultiply` | 80 × 200 and 80 × 180 dense | 0.01 ms | 0.00 ms | **1.52×** |
| `CountAnchorSharedNeighbors` | 80 anchors, k = 10 | 0.15 ms | 0.01 ms | **13.77×** |
| `EigenGramPCA` | 120 features × 400 cells, 10 PCs | 0.92 ms | 0.72 ms | **1.28×** |
| `SCTResidualStatsAndCorrected` | 800 × 400 sparse SCT | 2.34 ms | 2.25 ms | **1.04×** |
| `SCTPearsonResidualMatrix` | 800 × 400 sparse, 80 features | 0.41 ms | 0.15 ms | **2.78×** |
| `fast_dist` | 800 × 20 dense, 10 neighbors | 0.22 ms | 0.17 ms | **1.29×** |
| `row_sum_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.12 ms | 0.13 ms | 0.92× |
| `row_mean_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.12 ms | 0.12 ms | 0.98× |
| `row_var_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.46 ms | 0.41 ms | **1.10×** |

Times are medians from [`microbenchmark`](https://cran.r-project.org/package=microbenchmark) after one untimed warmup (3 to 11 repeats; fewer repeats on the slower kernels). Inputs are built before the timer starts. Reproduce the full table with Seurat and seurust installed:

```sh
Rscript scripts/bench-all-kernels.R
```

### Numeric error

Same inputs as the timing table. Each row is the difference between the Seurat C++ result and the seurust result.

Relative error is `|a − b| / max(|a|, |b|)` on entries with magnitude at least 1e-8. Significant figures are `−log10` of the largest of those relative errors. **exact** means the results matched bit for bit. Sparse matrices are compared by coordinate, and a missing entry counts as zero. `RunUMISampling` and `RunUMISamplingPerCell` draw random samples; both sides were started with `set.seed(1)`. `EigenGramPCA` loadings/embeddings are compared after aligning column signs (Spectra vs nalgebra can differ by a global sign per PC).

Most kernels match exactly or to 11+ figures. `EigenGramPCA` is the widest gap (7.2 figures): Seurat uses Spectra's partial eigensolver; seurust uses a full symmetric eigen-decomposition. `SCTPearsonResidualMatrix` agrees to 9.2 figures.

| Function | Max absolute error | Max relative error | Significant figures |
| --- | ---: | ---: | ---: |
| `LogNorm` | 8.88e-16 | 2.21e-16 | **15.7** |
| `Standardize` | 8.88e-16 | 2.83e-12 | **11.5** |
| `FastCov` | 8.88e-16 | 1.62e-12 | **11.8** |
| `FastCovMats` | 1.04e-16 | 3.23e-13 | **12.5** |
| `FastRBind` | 0 | 0 | exact |
| `RowVar` | 1.33e-15 | 1.10e-15 | **15.0** |
| `FastExpMean` | 2.66e-15 | 5.14e-16 | **15.3** |
| `SparseRowVar` | 3.55e-15 | 2.60e-15 | **14.6** |
| `SparseRowVarStd` | 2.66e-15 | 7.01e-16 | **15.2** |
| `SparseRowMeanVar` | 0 | 0 | exact |
| `SparseRowMeanVarLegacy` | 3.55e-15 | 2.60e-15 | **14.6** |
| `SparseRowVarStdLegacyRows` | 6.22e-15 | 1.47e-15 | **14.8** |
| `FastLogVMR` | 7.11e-15 | 5.86e-16 | **15.2** |
| `FastSparseRowScale` | 0 | 0 | exact |
| `FastDenseRowScale` | 0 | 0 | exact |
| `FastSparseRowScaleWithKnownStats` | 1.78e-15 | 2.22e-16 | **15.7** |
| `RowMergeMatrices` | 0 | 0 | exact |
| `ReplaceColsC` | 0 | 0 | exact |
| `GraphToNeighborHelper` | 0 | 0 | exact |
| `RunUMISampling` | 0 | 0 | exact |
| `RunUMISamplingPerCell` | 0 | 0 | exact |
| `ComputeSNN` | 0 | 0 | exact |
| `IntegrateDataC` | 0 | 0 | exact |
| `FindWeightsC` | 0 | 0 | exact |
| `ScoreHelper` | 2.22e-16 | 4.72e-16 | **15.3** |
| `SNN_SmallestNonzero_Dist` | 0 | 0 | exact |
| `RunModularityClusteringCpp` | 0 | 0 | exact |
| `RunModularityClusteringCpp_multi` | 0 | 0 | exact |
| `FindAllMarkersSparseFoldChangeStats` | 3.73e-09 | 7.00e-16 | **15.2** |
| `CcaCrossprodMultiply` | 0 | 0 | exact |
| `CountAnchorSharedNeighbors` | 0 | 0 | exact |
| `EigenGramPCA` | 8.34e-10 | 6.77e-08 | **7.2** |
| `SCTResidualStatsAndCorrected` | 3.08e-14 | 2.65e-11 | **10.6** |
| `SCTPearsonResidualMatrix` | 3.55e-15 | 6.04e-10 | **9.2** |
| `fast_dist` | 0 | 0 | exact |
| `row_sum_dgcmatrix` | 0 | 0 | exact |
| `row_mean_dgcmatrix` | 0 | 0 | exact |
| `row_var_dgcmatrix` | 0 | 0 | exact |

Reproduce with:

```sh
SEURUST_BENCH_MODE=error Rscript scripts/bench-all-kernels.R
```

### Machine

These numbers were measured on the machine that ran the script:

| | |
| --- | --- |
| OS | Ubuntu 24.04.5 LTS (x86_64), kernel 6.18.44 |
| Host | Cloud VM (shared virtualized host, not a dedicated benchmark machine) |
| CPU | Intel Xeon Processor @ 2.10 GHz, 4 vCPUs (1 thread per core); seurust and Seurat run with `setThreads(3)` |
| Memory | 15 GiB |
| R | 4.3.3 |
| Seurat | 5.6.0 |
| SeuratObject | 5.2.0 |
| Matrix | 1.6.5 |
| microbenchmark | 1.5.0 |
| seurust | 0.3.0 + unreleased #45 changes, commit `f1443e1` |
| Rust | rustc 1.97.0, release profile |

Shared cloud VMs add noise: run-to-run variation of roughly ±10% is expected, so ratios near 1.00× should be read as parity.

### End-to-end scRNA-seq workflow

This is a separate, earlier measurement on Ubuntu 22.04 (Docker dev image, R 4.6, Rust 1.95), not the machine in the table above, and it predates Seurat 5.6.0. It ran a simulated PBMC-style pipeline (~2,500 cells, 2,000 genes) with identical steps for both backends. PCA and UMAP used the same R code; only native kernel calls differed. **Speedup = C++ time ÷ Rust time** (> 1.0 means Rust is faster). The helper scripts for that workflow are not in this repository anymore. Kernel timings above are the ones `scripts/bench-all-kernels.R` still reproduces.

| Step | C++ (s) | Rust (s) | Speedup |
|------|--------:|---------:|--------:|
| QC native stats | 0.136 | 0.112 | **1.21×** |
| Log normalize | 0.054 | 0.054 | 1.00× |
| Variable features | 0.068 | 0.100 | 0.68× |
| Scale HVGs | 0.601 | 0.554 | **1.08×** |
| SNN graph | 0.992 | 0.942 | **1.05×** |
| Clustering | 0.234 | 0.280 | 0.84× |
| Batch integration | 0.309 | 0.284 | **1.09×** |
| **Total native kernels** | **2.39** | **2.33** | **1.03×** |

Cluster assignments, normalization digests, SNN structure, and integration outputs **matched exactly** between backends on that run. Summed native-kernel time was ~**3% faster** with Rust on this workflow (1.03× overall).

---

## What's ported

| Module | Seurat (upstream C++) | seurust | Status |
|--------|----------------------|---------|--------|
| Sparse row stats | `stats.cpp` | [`seurust/src/rust/src/stats.rs`](seurust/src/rust/src/stats.rs) | Ported |
| Data manipulation | `data_manipulation.cpp` | [`seurust/src/rust/src/data_manipulation/`](seurust/src/rust/src/data_manipulation/) | Ported (incl. 5.6 VST / scale / markers) |
| Integration | `integration.cpp` | [`seurust/src/rust/src/integration.rs`](seurust/src/rust/src/integration.rs) | Ported (incl. 5.6 CCA matvec) |
| SNN / kNN | `snn.cpp`, `fast_NN_dist.cpp` | [`snn.rs`](seurust/src/rust/src/snn.rs), [`fast_nn_dist.rs`](seurust/src/rust/src/fast_nn_dist.rs) | Ported |
| PCA | `pca.cpp` | [`seurust/src/rust/src/pca.rs`](seurust/src/rust/src/pca.rs) | Ported |
| SCTransform | `sctransform.cpp` | [`seurust/src/rust/src/sct.rs`](seurust/src/rust/src/sct.rs) | Ported |
| Modularity | `ModularityOptimizer.cpp` | [`seurust/src/rust/src/modularity.rs`](seurust/src/rust/src/modularity.rs) | Ported (incl. multi-resolution) |
| Annoy search | `annoy_nn.cpp` | — | Not ported |

The modularity optimizer is a Rust port of Seurat's C++ translation of Ludo Waltman's Java [ModularityOptimizer](http://www.ludowaltman.nl/slm/). A copy of that C++ translation is kept at [`seurust/src/cpp/`](seurust/src/cpp/) for reference. The package build links the Rust static library only.

Seurat 5.6.0 also added `AnnoySearchCpp` (on-disk Annoy queries). That entry point is not ported; `FindNeighbors` still uses Seurat's Annoy backend.

---

## Quick start

### Install seurust

Requires R ≥ 4.0 and a [Rust toolchain](https://rustup.rs) (rustc + Cargo ≥ 1.81) when installing from source.

**r-universe** or **GitHub**:

```r
install.packages(
  "seurust",
  repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
)

# Or development:
# remotes::install_github("NebilI/seurust", subdir = "seurust")
```

See [`seurust/README.md`](seurust/README.md) and [`r-universe/README.md`](r-universe/README.md) for packaging details.

### Verify parity in one line

```r
library(Seurat)
library(seurust)
library(Matrix)

mat <- Matrix::sparseMatrix(
  i = c(0, 2, 1), p = c(0, 1, 2, 3), x = 1:3, dims = c(3, 3)
)

mat_seurat <- mat
mat_seurat@x <- Seurat:::LogNorm(mat@x, mat@p, 1e4, 1L, FALSE)
all.equal(mat_seurat, seurust::LogNorm(mat, 1e4, FALSE))
# [1] TRUE
```

---

## Development

```sh
docker compose -f docker/docker-compose.yml build
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  bash docker/scripts/build-and-test-rust.sh
```

Full developer docs: [`docker/README.md`](docker/README.md).

Publishing the Rust crate and the r-universe R package is `release-seurust` (a push to `main` that changes `seurust/DESCRIPTION`, or **Actions → Release seurust**). That workflow creates the GitHub Release with notes from [`seurust/NEWS.md`](seurust/NEWS.md), publishes to crates.io, and updates `NebilI/NebilI.r-universe.dev`. How to cut a release: [`RELEASING.md`](RELEASING.md). Setup notes: [`r-universe/README.md`](r-universe/README.md).

---

## Relationship to Seurat

seurust follows native interfaces from [satijalab/seurat](https://github.com/satijalab/seurat). Upstream Seurat documentation still applies for analysis workflows:

- https://satijalab.org/seurat
- https://github.com/satijalab/seurat

`seurust/DESCRIPTION` lists Rahul Satija as a contributor for those original native APIs, and the Satija Lab and collaborators as copyright holders of the original Seurat / ModularityOptimizer code.

Contributions welcome — especially kernel ports, parity tests, and benchmark improvements. Open an [issue](https://github.com/NebilI/seurust/issues) or PR on this repository.

---

## License

MIT, the same license as current Seurat. See [LICENSE](LICENSE) and [LICENSE.md](LICENSE.md).

Seurat 5.6.0, the version these kernels were checked against, is `License: MIT + file LICENSE`. Seurat 4.0.0 was `GPL-3 | file LICENSE`; by Seurat 4.4.0 the upstream package had switched to MIT. This package stays on that MIT grant. Copyright for the original Seurat interfaces and the ModularityOptimizer translation is recorded in [`seurust/DESCRIPTION`](seurust/DESCRIPTION) and [`seurust/LICENSE`](seurust/LICENSE).

The modularity optimizer is a port of Seurat's C++ translation of Ludo Waltman's Java [ModularityOptimizer](http://www.ludowaltman.nl/slm/). That translation is included in MIT-licensed Seurat and the copy in this tree does not add a separate GPL header.
