# Docker development images

There is no Dockerfile in this repo for the upstream C++ tree; CI uses the pre-built [`satijalab/seurat-ci`](https://hub.docker.com/r/satijalab/seurat-ci) image. These images are for **local development** with your checkout mounted at `/workspace`.

| Image | Dockerfile | Purpose |
|-------|------------|---------|
| `rust-seurat-rcpp:dev` | `Dockerfile.rcpp` | R 4.x (r2u), compilers, Rcpp/RcppEigen/RcppProgress, SeuratObject |
| `rust-seurat-rust:dev` | `Dockerfile.rust` | Above + Rust stable, bindgen/clang, `rextendr` |

## Quick start

From the repository root:

```sh
docker compose -f docker/docker-compose.yml build
docker compose -f docker/docker-compose.yml run --rm rcpp-dev
```

Rust / extendr work:

```sh
docker compose -f docker/docker-compose.yml run --rm rust-dev
```

## Common tasks inside the container

**Compile and load the package (Rcpp + Rust path):**

```sh
Rscript -e "pkgbuild::compile_dll(debug = FALSE)"
Rscript -e "devtools::load_all()"
```

The Rust crate lives in `src/rust/`. `configure` generates `src/Makevars` from `src/Makevars.in` and runs `cargo build` + the `document` binary to refresh `R/extendr-wrappers.R`.

**Run Rust unit tests for the ported code:**

```sh
cargo test --manifest-path src/rust/Cargo.toml
```

**Regenerate Rcpp exports after editing remaining `src/*.cpp`:**

```sh
Rscript -e "Rcpp::compileAttributes()"
Rscript tools/fix-rcpp-init.R
Rscript -e "devtools::document()"
```

After `Rcpp::compileAttributes()`, always run `tools/fix-rcpp-init.R` so `entrypoint.c` remains the unified `R_init_Seurat` entry point for both Rcpp and extendr.

**Regenerate extendr wrappers after editing `src/rust/`:**

```sh
Rscript -e "rextendr::document()"
```

**Full dependency install** (heavy; matches more of CI):

```sh
Rscript /usr/local/bin/install-r-deps.R --full
```

**Standalone modularity optimizer** (optional sanity check of `ModularityOptimizer.cpp`):

```sh
cd src
clang++ -O3 -std=c++11 -DSTANDALONE -Wall -g ModularityOptimizer.cpp -o modularity_optimizer
```

**Check Rust / extendr setup:**

```sh
Rscript -e "rextendr::rust_sitrep()"
```

**End-to-end build + smoke test (Rust row stats):**

```sh
bash docker/scripts/build-and-test-rust.sh
```

**C++ vs Rust timing benchmarks** (informational; set `SEURAT_REQUIRE_RUST_FASTER=1` to fail when Rust is slower):

```sh
Rscript docker/scripts/benchmark-rust-cpp.R
```

Ratio `> 1.0` means Rust is faster. Modularity currently calls the same C++ optimizer through a bridge, so C++ is expected to win until a pure Rust port lands.

## Rust rewrite status

C++ remains the default for all production R code paths until you opt in to Rust.
Rust ports live alongside C++ with a `_rust` suffix on exported R functions.

| Module | C++ source (active) | Rust source (parallel) | Default R path |
|--------|---------------------|------------------------|----------------|
| Sparse row stats | `src/stats.cpp` | `src/rust/src/stats.rs` | `R/RcppExports.R` → `R/utilities.R` |
| Data manipulation | `src/data_manipulation.cpp` | `src/rust/src/data_manipulation/` | `R/RcppExports.R` |
| Integration | `src/integration.cpp` | `src/rust/src/integration.rs` | `R/RcppExports.R` |
| SNN / kNN | `src/snn.cpp`, `src/fast_NN_dist.cpp` | `src/rust/src/snn.rs`, `fast_nn_dist.rs` | `R/RcppExports.R` |
| Modularity | `src/ModularityOptimizer.cpp` | `src/rust/` (C++ bridge) | `R/RcppExports.R` |

Rust data-manipulation ports use **`ndarray`** (dense linear algebra) and **`sprs`** (sparse CSC/CSR).
Test helpers in `R/rust-sparse.R` wrap slot extraction/reconstruction for parity checks.

Compare C++ vs Rust for log-normalization:

```r
mat <- Matrix::sparseMatrix(i = c(0, 2, 1), p = c(0, 1, 2, 3), x = 1:3, dims = c(3, 3))
all.equal(LogNorm(mat, 1e4, FALSE), LogNormRust(mat, 1e4, FALSE))
```

## Build without Compose

```sh
docker build -f docker/Dockerfile.rcpp -t rust-seurat-rcpp:dev .
docker build -f docker/Dockerfile.rust -t rust-seurat-rust:dev .
docker run --rm -it -v "$(pwd):/workspace" -w /workspace rust-seurat-rcpp:dev
```

On Windows PowerShell, replace `$(pwd)` with `${PWD}`.

## Notes

- Base image `rocker/r2u:jammy` matches the **r2u** stack referenced in `.github/workflows/merge_checks.yaml` (Ubuntu + fast binary R packages).
- The rust image mounts a named volume at `src/rust/target` so `cargo` artifacts stay off the bind mount and out of git.
- Production/user-facing images remain [`satijalab/seurat`](https://hub.docker.com/r/satijalab/seurat); these are dev-only.
