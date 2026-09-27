# seurust

**A faster native backend for [Seurat](https://satijalab.org/seurat) — same workflows, same outputs, less time in the hot path.**

This repository is a development fork of Seurat v5 that adds **[seurust](seurust/)**, a companion R package with Rust/extendr reimplementations of Seurat's performance-critical native routines. Install both packages, keep your existing analysis code, and swap the backend where kernels have been ported.

> **Drop-in by design.** seurust exposes the same function signatures as Seurat's internal C++ layer (`LogNorm`, `FastSparseRowScale`, `ComputeSNN`, `IntegrateDataC`, and more). Parity tests assert bit-for-bit agreement with the original implementation on every ported routine.

[![seurust CI](https://github.com/NebilI/seurust/actions/workflows/seurust_checks.yaml/badge.svg)](https://github.com/NebilI/seurust/actions/workflows/seurust_checks.yaml)
[![r-universe](https://NebilI.r-universe.dev/badges/seurust)](https://NebilI.r-universe.dev/seurust)

Publishing notes: [`seurust/CRAN.md`](seurust/CRAN.md).
PR tests: `seurust Checks`. CRAN updates: Actions → **Build / submit seurust to CRAN**.

---

## Why use seurust?

Single-cell pipelines spend a surprising amount of time in a handful of native kernels: log-normalization, variable-feature statistics, scaling sparse matrices, building shared-nearest-neighbor graphs, and batch-integration weighting. Those routines dominate preprocessing and graph construction on large datasets.

seurust targets exactly that layer:

- **Same Seurat API** — no new object model, no workflow rewrite
- **Validated parity** — automated C++ vs Rust tests on every ported kernel ([`tests/testthat/test_rust_cpp_*.R`](tests/testthat/))
- **Measurable speedups** — Rust wins on the preprocessing kernels that run on every dataset (see benchmarks below)
- **Open development** — install from GitHub, run benchmarks locally, contribute kernel ports

Seurat itself remains the user-facing package. **seurust** is the engine upgrade you install alongside it.

---

## Architecture

```text
┌─────────────────────────────────────────────────────────────┐
│  Your R workflow (CreateSeuratObject, RunPCA, RunUMAP, …)   │
└───────────────────────────┬─────────────────────────────────┘
                            │
              ┌─────────────┴─────────────┐
              ▼                           ▼
     ┌─────────────────┐         ┌─────────────────┐
     │  Seurat (root)  │         │     seurust     │
     │  C++ / Rcpp     │         │  Rust / extendr │
     │  production API │         │  ported kernels │
     └─────────────────┘         └─────────────────┘
```

| Package | Location | Backend | Rust required? |
|---------|----------|---------|----------------|
| **Seurat** | repo root | C++/Rcpp | No |
| **seurust** | [`seurust/`](seurust/) | Rust/extendr | Yes (build time) |

---

## Performance vs Seurat (C++)

Median runtime for every ported native function, Seurat's C++ kernels against seurust. **Speedup = Seurat time ÷ seurust time.** Values **above 1.0× mean seurust is faster**. On the baseline run below, **29 of 30** kernels were faster; `ScoreHelper` was within timing noise (0.82× on this machine).

### seurust 0.1.0 (all 30 kernels)

Measured after syncing this fork to **Seurat 5.5.1** (see [Machine](#machine)).

| Function | Problem size | Seurat C++ | seurust Rust | Speedup |
| --- | --- | ---: | ---: | ---: |
| `LogNorm` | 5,000 × 2,000 sparse, 5% nonzero | 4.70 ms | 4.37 ms | **1.08×** |
| `Standardize` | 2,000 × 120 dense | 0.27 ms | 0.24 ms | **1.15×** |
| `FastCov` | 2,000 × 80 dense | 0.83 ms | 0.24 ms | **3.46×** |
| `FastCovMats` | 2,000 × 80 and 2,000 × 40 dense | 0.55 ms | 0.32 ms | **1.72×** |
| `FastRBind` | two 4,000 × 150 dense | 2.63 ms | 1.82 ms | **1.44×** |
| `RowVar` | 2,000 × 120 dense | 0.09 ms | 0.04 ms | **2.11×** |
| `FastExpMean` | 5,000 × 2,000 sparse, 5% nonzero | 4.10 ms | 0.62 ms | **6.62×** |
| `SparseRowVar` | 5,000 × 2,000 sparse, 5% nonzero | 2.47 ms | 0.40 ms | **6.12×** |
| `SparseRowVar2` | 5,000 × 2,000 sparse, 5% nonzero | 2.72 ms | 0.28 ms | **9.83×** |
| `SparseRowVarStd` | 5,000 × 2,000 sparse, 5% nonzero | 2.81 ms | 0.33 ms | **8.59×** |
| `FastLogVMR` | 5,000 × 2,000 sparse, 5% nonzero | 5.97 ms | 1.29 ms | **4.64×** |
| `FastSparseRowScale` | 5,000 × 2,000 sparse, 5% nonzero | 84.6 ms | 32.6 ms | **2.60×** |
| `FastSparseRowScaleWithKnownStats` | 5,000 × 2,000 sparse, 5% nonzero | 89.9 ms | 31.9 ms | **2.81×** |
| `RowMergeMatrices` | 1,500 × 400 and 1,200 × 350 CSR | 0.58 ms | 0.21 ms | **2.71×** |
| `ReplaceColsC` | 2,000 × 600 sparse, replace 20 columns | 0.31 ms | 0.09 ms | **3.55×** |
| `GraphToNeighborHelper` | 2,000 cells, 30 neighbors, symmetric | 0.52 ms | 0.44 ms | **1.18×** |
| `RunUMISampling` | 2,000 × 1,500 counts, 8% nonzero | 3.12 ms | 2.54 ms | **1.23×** |
| `RunUMISamplingPerCell` | 2,000 × 1,500 counts, 8% nonzero | 2.96 ms | 2.35 ms | **1.26×** |
| `ComputeSNN` | 1,500 cells, k = 20 | 4.41 ms | 2.09 ms | **2.11×** |
| `IntegrateDataC` | 800 genes × 200 cells, 60 anchors | 0.76 ms | 0.54 ms | **1.42×** |
| `FindWeightsC` | 400 cells, 600 anchors, k = 10 | 0.58 ms | 0.49 ms | **1.20×** |
| `ScoreHelper` | 300 cells × 16 dimensions | 0.71 ms | 0.87 ms | 0.82× |
| `WriteEdgeFile` | 400-cell SNN | 1.60 ms | 1.43 ms | **1.12×** |
| `DirectSNNToFile` | 400 cells, k = 15 | 2.49 ms | 1.84 ms | **1.35×** |
| `SNN_SmallestNonzero_Dist` | 400 cells × 10 dimensions | 0.43 ms | 0.33 ms | **1.31×** |
| `RunModularityClusteringCpp` | 250-cell SNN, 10 iterations | 3.96 ms | 3.81 ms | **1.04×** |
| `fast_dist` | 800 × 20 dense, 10 neighbors | 0.22 ms | 0.20 ms | **1.13×** |
| `row_sum_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.21 ms | 0.11 ms | **2.04×** |
| `row_mean_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.23 ms | 0.11 ms | **2.14×** |
| `row_var_dgcmatrix` | 6,000 × 1,500 sparse, 4% nonzero | 0.74 ms | 0.34 ms | **2.16×** |

Times are medians from [`microbenchmark`](https://cran.r-project.org/package=microbenchmark) after one untimed warmup (3 to 11 repeats; fewer repeats on the slower kernels). Inputs are built before the timer starts. Reproduce the full table with:

```sh
Rscript scripts/bench-all-kernels.R
```

### Numeric error

Same inputs as the timing table. Each row is the difference between the Seurat C++ result and the seurust result.

Relative error is `|a − b| / max(|a|, |b|)` on entries with magnitude at least 1e-8. Significant figures are `−log10` of the largest of those relative errors. **exact** means the results matched bit for bit. Sparse matrices are compared by coordinate, and a missing entry counts as zero. `WriteEdgeFile` is compared on the edge file it writes. `DirectSNNToFile` is compared on both the returned matrix and the edge file. `RunUMISampling` and `RunUMISamplingPerCell` draw random samples; both sides were started with `set.seed(1)`.

18 of the 30 kernels match exactly. The other twelve agree to at least 11.5 significant figures. `Standardize` has the widest gap (relative error 2.83e-12). `FastCov` and `FastCovMats` agree to about 12 figures. `FastExpMean`, `FastLogVMR`, `SparseRowVar`, `SparseRowVarStd`, `SparseRowVar2`, `RowVar`, and `ScoreHelper` agree to 14.7–15.3 figures. The two sparse row-scaling kernels differ by one unit in the last place (relative error 2.22e-16). `FastExpMean`, `FastLogVMR`, `SparseRowVar`, and `SparseRowVarStd` are no longer bit-identical to Seurat: the parallel column scan adds values in a different order.

| Function | Max absolute error | Max relative error | Significant figures |
| --- | ---: | ---: | ---: |
| `LogNorm` | 0 | 0 | exact |
| `Standardize` | 5.77e-15 | 2.83e-12 | **11.5** |
| `FastCov` | 8.88e-16 | 1.62e-12 | **11.8** |
| `FastCovMats` | 1.04e-16 | 3.23e-13 | **12.5** |
| `FastRBind` | 0 | 0 | exact |
| `RowVar` | 1.33e-15 | 1.10e-15 | **15.0** |
| `FastExpMean` | 2.66e-15 | 6.02e-16 | **15.2** |
| `SparseRowVar` | 3.11e-15 | 2.18e-15 | **14.7** |
| `SparseRowVar2` | 3.11e-15 | 2.18e-15 | **14.7** |
| `SparseRowVarStd` | 7.99e-15 | 1.94e-15 | **14.7** |
| `FastLogVMR` | 7.11e-15 | 5.86e-16 | **15.2** |
| `FastSparseRowScale` | 1.78e-15 | 2.22e-16 | **15.7** |
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
| `WriteEdgeFile` | 0 | 0 | exact |
| `DirectSNNToFile` | 0 | 0 | exact |
| `SNN_SmallestNonzero_Dist` | 0 | 0 | exact |
| `RunModularityClusteringCpp` | 0 | 0 | exact |
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
| OS | Ubuntu 24.04.4 LTS (Linux 6.12.94+, x86_64) |
| CPU | 4 cores |
| Memory | 15 GiB |
| R | 4.3.3 |
| Seurat | 5.5.1 (this fork, synced with upstream) |
| seurust | 0.1.0 |
| Rust | rustc 1.98.1 |

### End-to-end scRNA-seq workflow

This is a separate, earlier measurement on Ubuntu 22.04 (Docker dev image, R 4.6, Rust 1.95), not the machine in the table above. It runs a simulated PBMC-style pipeline (~2,500 cells, 2,000 genes) with identical steps for both backends ([`examples/compare_scrna_workflows.R`](examples/compare_scrna_workflows.R)). PCA and UMAP use the same R code; only native kernel calls differ. **Speedup = C++ time ÷ Rust time** (> 1.0 means Rust is faster).

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

Cluster assignments, normalization digests, SNN structure, and integration outputs **match exactly** between backends. Summed native-kernel time is ~**3% faster** with Rust on this workflow (1.03× overall).

Reproduce locally:

```sh
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  Rscript docker/scripts/benchmark-rust-cpp.R

docker compose -f docker/docker-compose.yml run --rm rust-dev \
  Rscript examples/compare_scrna_workflows.R
```

---

## What's ported

| Module | Seurat (C++) | seurust (Rust) | Status |
|--------|--------------|----------------|--------|
| Sparse row stats | `src/stats.cpp` | `stats.rs` | ✅ Ported |
| Data manipulation | `src/data_manipulation.cpp` | `data_manipulation/` | ✅ Ported |
| Integration | `src/integration.cpp` | `integration.rs` | ✅ Ported |
| SNN / kNN | `src/snn.cpp`, `fast_NN_dist.cpp` | `snn.rs`, `fast_nn_dist.rs` | ✅ Ported |
| Modularity | `src/ModularityOptimizer.cpp` | C++ bridge | 🔶 Bridge (pure Rust port planned) |

---

## Quick start

### Install seurust

Requires R ≥ 4.0 and a [Rust toolchain](https://rustup.rs) (rustc + Cargo ≥ 1.81) when installing from source.

**CRAN** (after acceptance) or **r-universe / GitHub**:

```r
install.packages("seurust")

# Until CRAN accepts:
# install.packages(
#   "seurust",
#   repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
# )

# Or development:
# remotes::install_github("NebilI/seurust", subdir = "seurust")
```

See [`seurust/README.md`](seurust/README.md) and [`seurust/CRAN.md`](seurust/CRAN.md) for packaging details.

### Verify parity in one line

```r
library(Seurat)
library(seurust)
library(Matrix)

mat <- Matrix::sparseMatrix(
  i = c(0, 2, 1), p = c(0, 1, 2, 3), x = 1:3, dims = c(3, 3)
)

all.equal(
  Seurat:::LogNorm(mat, 1e4, FALSE),
  seurust::LogNorm(mat, 1e4, FALSE)
)
# [1] TRUE
```

### Run the example workflow

```sh
Rscript examples/scrna_workflow_rust.R    # Rust backend
Rscript examples/scrna_workflow_cpp.R     # C++ backend (comparison)
Rscript examples/compare_scrna_workflows.R
```

---

## Development

```sh
# Build dev environment (Docker recommended)
docker compose -f docker/docker-compose.yml build
docker compose -f docker/docker-compose.yml run --rm rust-dev

# Inside the container: compile, test, benchmark
Rscript docker/scripts/build-and-test-rust.sh
Rscript docker/scripts/benchmark-rust-cpp.R
```

Full developer docs: [`docker/README.md`](docker/README.md).

---

## Relationship to upstream Seurat

This fork tracks [satijalab/seurat](https://github.com/satijalab/seurat) and adds the Rust migration layer under [`seurust/`](seurust/). Upstream Seurat documentation and vignettes still apply for analysis workflows:

- https://satijalab.org/seurat
- https://cran.r-project.org/package=Seurat

Contributions welcome — especially kernel ports, parity tests, and benchmark improvements. Open an [issue](https://github.com/NebilI/seurust/issues) or PR on this repository.

---

## License

MIT — same as upstream Seurat. See [LICENSE](LICENSE).
