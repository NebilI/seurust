# seurust

Rust/extendr backend for Seurat's performance-critical native routines. Install
alongside [Seurat](../) to compare C++ and Rust implementations during the
migration.

## Requirements

- R (>= 4.0.0)
- Rust toolchain: [rustc](https://www.rust-lang.org/tools/install) and Cargo (>= 1.65)
- On Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) plus Rust

## Install

### From GitHub (recommended)

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("NebilI/seurust", subdir = "seurust")
```

Install a specific branch or tag:

```r
remotes::install_github(
  "NebilI/seurust",
  subdir = "seurust",
  ref = "feature/rust-rewrite"
)
```

### From a release tarball

Download `seurust_*.tar.gz` from [GitHub Releases](https://github.com/NebilI/seurust/releases), then:

```r
install.packages("path/to/seurust_0.1.0.tar.gz", repos = NULL, type = "source")
```

### From r-universe

After registering this package in your [r-universe](https://r-universe.dev) registry (see [`r-universe/`](../r-universe/README.md)):

```r
install.packages(
  "seurust",
  repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
)
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

mat <- Matrix::sparseMatrix(i = c(0, 2, 1), p = c(0, 1, 2, 3), x = 1:3, dims = c(3, 3))
all.equal(
  Seurat:::LogNorm(mat, 1e4, FALSE),
  seurust::LogNorm(mat, 1e4, FALSE)
)
```

Parity and benchmark tests live in the parent package under
`tests/testthat/test_rust_cpp_*.R` and require `seurust` in `Suggests`.

## Layout

| Path | Role |
|------|------|
| `src/rust/` | extendr crate (Rust kernels) |
| `src/cpp/` | ModularityOptimizer C++ bridge |
| `src/entrypoint.c` | Links Rust staticlib into `seurust.so` |
| `R/native.R` | High-level R API matching Seurat's RcppExports |
| `R/extendr-wrappers.R` | Generated low-level `.Call` wrappers |

Seurat itself is C++/Rcpp-only; no Rust toolchain is required to build the main
package.
