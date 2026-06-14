# Rust-Seurat

**A faster native backend for [Seurat](https://satijalab.org/seurat) — same workflows, same outputs, less time in the hot path.**

This repository is a development fork of Seurat v5 that adds **[RSeurat](RSeurat/)**, a companion R package with Rust/extendr reimplementations of Seurat's performance-critical native routines. Install both packages, keep your existing analysis code, and swap the backend where kernels have been ported.

> **Drop-in by design.** RSeurat exposes the same function signatures as Seurat's internal C++ layer (`LogNorm`, `FastSparseRowScale`, `ComputeSNN`, `IntegrateDataC`, and more). Parity tests assert bit-for-bit agreement with the original implementation on every ported routine.

[![RSeurat CI](https://github.com/NebilI/Rust-Seurat/actions/workflows/rseurat_checks.yaml/badge.svg)](https://github.com/NebilI/Rust-Seurat/actions/workflows/rseurat_checks.yaml)

---

## Why use Rust-Seurat?

Single-cell pipelines spend a surprising amount of time in a handful of native kernels: log-normalization, variable-feature statistics, scaling sparse matrices, building shared-nearest-neighbor graphs, and batch-integration weighting. Those routines dominate preprocessing and graph construction on large datasets.

Rust-Seurat targets exactly that layer:

- **Same Seurat API** — no new object model, no workflow rewrite
- **Validated parity** — automated C++ vs Rust tests on every ported kernel ([`tests/testthat/test_rust_cpp_*.R`](tests/testthat/))
- **Measurable speedups** — Rust wins on the preprocessing kernels that run on every dataset (see benchmarks below)
- **Open development** — install from GitHub, run benchmarks locally, contribute kernel ports

Seurat itself remains the user-facing package. **RSeurat** is the engine upgrade you install alongside it.

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
     │  Seurat (root)  │         │     RSeurat     │
     │  C++ / Rcpp     │         │  Rust / extendr │
     │  production API │         │  ported kernels │
     └─────────────────┘         └─────────────────┘
```

| Package | Location | Backend | Rust required? |
|---------|----------|---------|----------------|
| **Seurat** | repo root | C++/Rcpp | No |
| **RSeurat** | [`RSeurat/`](RSeurat/) | Rust/extendr | Yes (build time) |

---

## Performance vs Seurat (C++)

Benchmarks below were collected on **Ubuntu 22.04** (Docker dev image, R 4.6, Rust 1.95) using [`docker/scripts/benchmark-rust-cpp.R`](docker/scripts/benchmark-rust-cpp.R). Each row is the median-time ratio **C++ ÷ Rust**; values **> 1.0 mean Rust is faster**.

### Kernel micro-benchmarks

| Routine | Problem size | Rust vs C++ | Winner |
|---------|--------------|-------------|--------|
| **FastSparseRowScale** | 2,000 × 2,500 sparse | **1.40×** | Rust |
| **LogNorm** | 400 × 400 sparse | **1.33×** | Rust |
| **SparseRowVar2** | 2,000 × 2,500 sparse | **1.23×** | Rust |
| **row_sum_dgcmatrix** | 3,000 × 800 sparse | **2.86×** | Rust |
| Modularity clustering | 34-node SNN, 5×50 iters | 1.00× | Tie |
| **ComputeSNN** | 500 cells, *k* = 20 | 0.60× | C++ |
| **ComputeSNN** | 2,000 cells, *k* = 20 | 0.95× | ~Tie |

**Takeaway:** Rust delivers the largest gains on sparse matrix preprocessing — normalization, scaling, and row statistics — the steps that run on every dataset. SNN graph construction is actively being optimized; at 2,000 cells the backends are already within ~5%.

### End-to-end scRNA-seq workflow

We run a full simulated PBMC-style pipeline (~2,500 cells, 2,000 genes) with identical steps for both backends ([`examples/compare_scrna_workflows.R`](examples/compare_scrna_workflows.R)). PCA and UMAP use the same R code; only native kernel calls differ. **Speedup = C++ time ÷ Rust time** (same convention as above; > 1.0 means Rust is faster).

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

| Module | Seurat (C++) | RSeurat (Rust) | Status |
|--------|--------------|----------------|--------|
| Sparse row stats | `src/stats.cpp` | `stats.rs` | ✅ Ported |
| Data manipulation | `src/data_manipulation.cpp` | `data_manipulation/` | ✅ Ported |
| Integration | `src/integration.cpp` | `integration.rs` | ✅ Ported |
| SNN / kNN | `src/snn.cpp`, `fast_NN_dist.cpp` | `snn.rs`, `fast_nn_dist.rs` | ✅ Ported |
| Modularity | `src/ModularityOptimizer.cpp` | C++ bridge | 🔶 Bridge (pure Rust port planned) |

---

## Quick start

### Install RSeurat

Requires R ≥ 4.0 and a [Rust toolchain](https://rustup.rs) (rustc + Cargo ≥ 1.65).

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Rust backend
remotes::install_github("NebilI/Rust-Seurat", subdir = "RSeurat")

# Seurat from this fork (or use CRAN Seurat for comparison)
remotes::install_github("NebilI/Rust-Seurat")
```

See [`RSeurat/README.md`](RSeurat/README.md) for local build instructions and [r-universe setup](r-universe/README.md) for publishing pre-built packages.

### Verify parity in one line

```r
library(Seurat)
library(RSeurat)
library(Matrix)

mat <- Matrix::sparseMatrix(
  i = c(0, 2, 1), p = c(0, 1, 2, 3), x = 1:3, dims = c(3, 3)
)

all.equal(
  Seurat:::LogNorm(mat, 1e4, FALSE),
  RSeurat::LogNorm(mat, 1e4, FALSE)
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

This fork tracks [satijalab/seurat](https://github.com/satijalab/seurat) and adds the Rust migration layer under [`RSeurat/`](RSeurat/). Upstream Seurat documentation and vignettes still apply for analysis workflows:

- https://satijalab.org/seurat
- https://cran.r-project.org/package=Seurat

Contributions welcome — especially kernel ports, parity tests, and benchmark improvements. Open an [issue](https://github.com/NebilI/Rust-Seurat/issues) or PR on this repository.

---

## License

MIT — same as upstream Seurat. See [LICENSE](LICENSE).
