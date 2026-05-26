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

**Compile and load the package (Rcpp path):**

```sh
Rscript -e "pkgbuild::compile_dll(debug = FALSE)"
Rscript -e "devtools::load_all()"
```

**Full dependency install** (heavy; matches more of CI):

```sh
Rscript /usr/local/bin/install-r-deps.R --full
```

**Regenerate Rcpp exports after editing `src/*.cpp`:**

```sh
Rscript -e "Rcpp::compileAttributes()"
Rscript -e "devtools::document()"
```

**Standalone modularity optimizer** (optional sanity check of `ModularityOptimizer.cpp`):

```sh
cd src
clang++ -O3 -std=c++11 -DSTANDALONE -Wall -g ModularityOptimizer.cpp -o modularity_optimizer
```

**Check Rust / extendr setup:**

```sh
rextendr::rust_sitrep()
```

**Scaffold extendr in the package** (when you start the rewrite):

```sh
Rscript -e "rextendr::use_extendr()"
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
- The rust image sets `CARGO_TARGET_DIR=/workspace/src/rust-target` so `cargo` artifacts stay on a named volume and do not pollute git-tracked paths.
- Production/user-facing images remain [`satijalab/seurat`](https://hub.docker.com/r/satijalab/seurat); these are dev-only.
