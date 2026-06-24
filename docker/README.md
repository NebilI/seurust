# Docker development images

There is no Dockerfile in this repo for the upstream C++ tree; CI uses the pre-built [`satijalab/seurat-ci`](https://hub.docker.com/r/satijalab/seurat-ci) image. These images are for **local development** with your checkout mounted at `/workspace`.

| Image | Dockerfile | Purpose |
|-------|------------|---------|
| `seurust-rcpp:dev` | `Dockerfile.rcpp` | R 4.x (r2u), compilers, Rcpp/RcppEigen/RcppProgress, SeuratObject |
| `seurust-rust:dev` | `Dockerfile.rust` | Above + Rust stable, bindgen/clang, `rextendr` |

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

## Two-package layout

| Package | Native backend | Rust required? |
|---------|----------------|----------------|
| **Seurat** (root) | C++/Rcpp | No |
| **seurust** (`seurust/`) | Rust/extendr | Yes |

Install both to compare implementations:

```r
devtools::load_all()                    # Seurat (C++)
devtools::install("seurust")       # Rust backend
library(seurust)

all.equal(
  Seurat:::LogNorm(mat, 1e4, FALSE),
  seurust::LogNorm(mat, 1e4, FALSE)
)
```

## Common tasks inside the container

**Compile Seurat (C++ only):**

```sh
Rscript -e "pkgbuild::compile_dll(debug = FALSE)"
Rscript -e "devtools::load_all()"
```

**Build and install seurust:**

```sh
cd seurust && Rscript tools/config.R && cd ..
R CMD INSTALL seurust
```

The Rust crate lives in `seurust/src/rust/`. `seurust/configure` generates `Makevars` and runs `cargo build` + the `document` binary to refresh `R/extendr-wrappers.R`.

**Run Rust unit tests:**

```sh
cargo test --manifest-path seurust/src/rust/Cargo.toml
```

**Regenerate Rcpp exports after editing `src/*.cpp`:**

```sh
Rscript -e "Rcpp::compileAttributes()"
```

**Regenerate extendr wrappers after editing `seurust/src/rust/`:**

```sh
cd seurust/src/rust && cargo run --bin document --release && cd ../../..
```

**End-to-end build + parity checks** (installs Seurat Imports, then builds both packages):

```sh
bash docker/scripts/build-and-test-rust.sh
```

Or from the host in one shot:

```sh
docker compose -f docker/docker-compose.yml run --rm rust-dev bash docker/scripts/build-and-test-rust.sh
```

The first run installs many R packages from `DESCRIPTION` and can take several minutes.

Each `docker compose run` starts a **fresh container**. Standalone scripts call
`docker/scripts/bootstrap-dev-env.R` first to install Imports, compile Seurat,
and install seurust if needed:

```sh
docker compose -f docker/docker-compose.yml run --rm rust-dev Rscript docker/scripts/run-rust-parity.R
```

**C++ vs Rust timing benchmarks** (set `SEURAT_REQUIRE_RUST_FASTER=1` to fail when Rust is slower):

```sh
Rscript docker/scripts/benchmark-rust-cpp.R
```

Ratio `> 1.0` means Rust is faster. Modularity currently calls the same C++ optimizer through a bridge, so C++ is expected to win until a pure Rust port lands.

## Rust rewrite status

Seurat remains the production package (C++/Rcpp). **seurust** is a sibling package with the same R API for ported kernels, used for parity testing and benchmarks.

| Module | Seurat (C++) | seurust |
|--------|--------------|------------|
| Sparse row stats | `src/stats.cpp` | `seurust/src/rust/src/stats.rs` |
| Data manipulation | `src/data_manipulation.cpp` | `seurust/src/rust/src/data_manipulation/` |
| Integration | `src/integration.cpp` | `seurust/src/rust/src/integration.rs` |
| SNN / kNN | `src/snn.cpp`, `src/fast_NN_dist.cpp` | `seurust/src/rust/src/snn.rs`, `fast_nn_dist.rs` |
| Modularity | `src/ModularityOptimizer.cpp` | C++ bridge in `seurust/src/rust/` |

## Build without Compose

```sh
docker build -f docker/Dockerfile.rcpp -t seurust-rcpp:dev .
docker build -f docker/Dockerfile.rust -t seurust-rust:dev .
docker run --rm -it -v "$(pwd):/workspace" -w /workspace seurust-rust:dev
```

On Windows PowerShell, replace `$(pwd)` with `${PWD}`.

## Notes

- Base image `rocker/r2u:jammy` matches the **r2u** stack referenced in `.github/workflows/merge_checks.yaml`.
- The rust image mounts a named volume at `seurust/src/rust/target` so `cargo` artifacts stay off the bind mount.
- Production/user-facing images remain [`satijalab/seurat`](https://hub.docker.com/r/satijalab/seurat); these are dev-only.
