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

### seurust 0.2.0 (Seurat 5.6.0 kernels)

Measured against **Seurat 5.6.0** with `setThreads(3)` on this 4-core machine (see [Machine](#machine)), after the Rust kernel work in [`OPTIMIZATION_ROADMAP.md`](OPTIMIZATION_ROADMAP.md) (items 1–7). Seurat 5.6.0 dropped `SparseRowVar2`, `WriteEdgeFile`, and `DirectSNNToFile` from its R exports (seurust still implements them, but they are not in this comparison). `AnnoySearchCpp` is not ported. **34 of 38** kernels were faster. The four that still lose are near-parity or sub-millisecond: `CcaCrossprodMultiply` (0.71×), `EigenGramPCA` (0.78×), `SparseRowMeanVar` (0.80×), and `fast_dist` (0.90×).

| Function | Problem size | Seurat C++ | seurust Rust | Speedup |
| --- | --- | ---: | ---: | ---: |
| `LogNorm` | 5,000 × 2,000 sparse, 5% nonzero | 1.04 ms | 0.78 ms | **1.33×** |
| `Standardize` | 2,000 × 120 dense | 0.18 ms | 0.06 ms | **2.78×** |
| `FastCov` | 2,000 × 80 dense | 0.87 ms | 0.23 ms | **3.72×** |
| `FastCovMats` | 2,000 × 80 and 2,000 × 40 dense | 0.52 ms | 0.37 ms | **1.41×** |
| `FastRBind` | two 4,000 × 150 dense | 2.48 ms | 1.85 ms | **1.34×** |
| `RowVar` | 2,000 × 120 dense | 0.09 ms | 0.04 ms | **2.11×** |
| `FastExpMean` | 5,000 × 2,000 sparse, 5% nonzero | 4.06 ms | 0.55 ms | **7.32×** |
| `SparseRowVar` | 5,000 × 2,000 sparse, 5% nonzero | 2.45 ms | 0.42 ms | **5.83×** |
| `SparseRowVarStd` | 5,000 × 2,000 sparse, 5% nonzero | 0.59 ms | 0.30 ms | **1.95×** |
| `SparseRowMeanVar` | 5,000 × 2,000 sparse, 5% nonzero | 0.21 ms | 0.26 ms | 0.80× |
| `SparseRowMeanVarLegacy` | 5,000 × 2,000 sparse, 5% nonzero | 0.64 ms | 0.37 ms | **1.75×** |
| `SparseRowVarStdLegacyRows` | 5,000 × 2,000 sparse, 80 rows | 0.41 ms | 0.10 ms | **4.07×** |
| `FastLogVMR` | 5,000 × 2,000 sparse, 5% nonzero | 5.92 ms | 1.14 ms | **5.19×** |
| `FastSparseRowScale` | 5,000 × 2,000 sparse, 5% nonzero | 10.9 ms | 5.14 ms | **2.12×** |
| `FastDenseRowScale` | 800 × 300 dense | 0.40 ms | 0.18 ms | **2.20×** |
| `FastSparseRowScaleWithKnownStats` | 5,000 × 2,000 sparse, 5% nonzero | 90.6 ms | 5.28 ms | **17.16×** |
| `RowMergeMatrices` | 1,500 × 400 and 1,200 × 350 CSR | 0.56 ms | 0.18 ms | **3.10×** |
| `ReplaceColsC` | 2,000 × 600 sparse, replace 20 columns | 0.30 ms | 0.08 ms | **3.60×** |
| `GraphToNeighborHelper` | 2,000 cells, 30 neighbors, symmetric | 0.51 ms | 0.38 ms | **1.34×** |
| `RunUMISampling` | 2,000 × 1,500 counts, 8% nonzero | 2.90 ms | 2.18 ms | **1.33×** |
| `RunUMISamplingPerCell` | 2,000 × 1,500 counts, 8% nonzero | 2.84 ms | 2.08 ms | **1.37×** |
| `ComputeSNN` | 1,500 cells, k = 20 | 0.99 ms | 0.69 ms | **1.44×** |
| `IntegrateDataC` | 800 genes × 200 cells, 60 anchors | 0.83 ms | 0.38 ms | **2.17×** |
| `FindWeightsC` | 400 cells, 600 anchors, k = 10 | 0.26 ms | 0.12 ms | **2.18×** |
| `ScoreHelper` | 300 cells × 16 dimensions | 0.42 ms | 0.40 ms | **1.04×** |
| `SNN_SmallestNonzero_Dist` | 400 cells × 10 dimensions | 0.36 ms | 0.32 ms | **1.11×** |
| `RunModularityClusteringCpp` | 250-cell SNN, 10 iterations | 3.80 ms | 3.39 ms | **1.12×** |
| `RunModularityClusteringCpp_multi` | 250-cell SNN, 3 resolutions | 1.97 ms | 1.06 ms | **1.85×** |
| `FindAllMarkersSparseFoldChangeStats` | 5,000 × 2,000 sparse, 4 groups | 1.28 ms | 1.00 ms | **1.29×** |
| `CcaCrossprodMultiply` | 80 × 200 and 80 × 180 dense | 0.01 ms | 0.01 ms | 0.71× |
| `CountAnchorSharedNeighbors` | 80 anchors, k = 10 | 0.12 ms | 0.01 ms | **12.32×** |
| `EigenGramPCA` | 120 features × 400 cells, 10 PCs | 0.63 ms | 0.81 ms | 0.78× |
| `SCTResidualStatsAndCorrected` | 800 × 400 sparse SCT | 1.85 ms | 1.28 ms | **1.45×** |
| `SCTPearsonResidualMatrix` | 800 × 400 sparse, 80 features | 0.42 ms | 0.14 ms | **3.06×** |
| `fast_dist` | 800 × 20 dense, 10 neighbors | 0.18 ms | 0.19 ms | 0.90× |
| `row_sum_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.21 ms | 0.10 ms | **2.10×** |
| `row_mean_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.23 ms | 0.10 ms | **2.22×** |
| `row_var_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.74 ms | 0.34 ms | **2.16×** |

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
| `FastExpMean` | 2.66e-15 | 6.02e-16 | **15.2** |
| `SparseRowVar` | 3.11e-15 | 2.18e-15 | **14.7** |
| `SparseRowVarStd` | 2.66e-15 | 6.46e-16 | **15.2** |
| `SparseRowMeanVar` | 0 | 0 | exact |
| `SparseRowMeanVarLegacy` | 3.11e-15 | 2.18e-15 | **14.7** |
| `SparseRowVarStdLegacyRows` | 4.44e-15 | 1.17e-15 | **14.9** |
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
| `FindAllMarkersSparseFoldChangeStats` | 3.73e-09 | 6.20e-16 | **15.2** |
| `CcaCrossprodMultiply` | 0 | 0 | exact |
| `CountAnchorSharedNeighbors` | 0 | 0 | exact |
| `EigenGramPCA` | 8.33e-10 | 6.77e-08 | **7.2** |
| `SCTResidualStatsAndCorrected` | 1.14e-14 | 5.58e-12 | **11.3** |
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
| OS | Ubuntu 24.04 (x86_64) |
| CPU | 4 cores (`setThreads(3)`) |
| Memory | 15 GiB |
| R | 4.6.1 |
| Seurat | 5.6.0 |
| seurust | 0.2.0 |
| Rust | rustc 1.83.0 |

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
