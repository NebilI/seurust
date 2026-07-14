# seurust

Rust/extendr backend for Seurat's performance-critical native routines. Install
alongside [Seurat](https://satijalab.org/seurat) to use the same function
signatures with a Rust backend.

## Requirements

- R (>= 4.0.0)
- Rust toolchain: [rustc](https://rust-lang.org/tools/install/) and Cargo (>= 1.81)
- On Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) plus Rust

## Install

### From r-universe / GitHub

CRAN publication is on hold until the maintainer explicitly approves a release.

```r
install.packages(
  "seurust",
  repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
)
```

### From GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("NebilI/seurust", subdir = "seurust")
```

### From a release tarball

Download `seurust_*.tar.gz` from [GitHub Releases](https://github.com/NebilI/seurust/releases), then:

```r
install.packages("path/to/seurust_0.1.0.tar.gz", repos = NULL, type = "source")
```

### Local development

From the repo root:

```r
devtools::install("seurust")
```

Or from the shell:

```sh
cd seurust
Rscript tools/config.R
cd ..
R CMD INSTALL seurust
```

## Compare against Seurat

```r
library(Seurat)
library(seurust)
library(Matrix)

mat <- Matrix::sparseMatrix(i = c(1, 3, 2), j = c(1, 2, 3), x = 1:3, dims = c(3, 3))
all.equal(
  Seurat:::LogNorm(mat, 1e4, FALSE),
  seurust::LogNorm(mat, 1e4, FALSE)
)
```

## Publishing

See **[`CRAN.md`](CRAN.md)** for the full r-universe / GitHub Release / CRAN process
(all local steps use Docker Compose).

```sh
docker compose -f docker/docker-compose.yml run --rm seurust-cran
docker compose -f docker/docker-compose.yml run --rm -e SUBMIT_CRAN=yes seurust-cran-submit
```

| Workflow | Channel |
|----------|---------|
| `publish-seurust-r.yaml` | GitHub Release tarball + r-universe registry sync |
| `publish-seurust-crate.yaml` | crates.io (`seurust` crate) |
| `build-seurust-cran.yaml` | Docker Compose CRAN tarball + `--as-cran` check |

## Layout

| Path | Role |
|------|------|
| `src/rust/` | extendr crate (Rust kernels) |
| `src/rust/vendor.tar.xz` | Vendored crates for offline/CRAN builds |
| `src/cpp/` | ModularityOptimizer C++ bridge |
| `src/entrypoint.c` | Links Rust staticlib into `seurust.so` |
| `R/native.R` | High-level R API matching Seurat's RcppExports |
| `R/extendr-wrappers.R` | Generated low-level `.Call` wrappers |
