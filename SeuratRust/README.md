# SeuratRust

Rust/extendr backend for Seurat's performance-critical native routines. Install
alongside [Seurat](../) to compare C++ and Rust implementations during the
migration.

## Install

From the repo root (requires Rust toolchain):

```r
devtools::install("SeuratRust")
```

Or from the shell:

```sh
cd SeuratRust
Rscript tools/config.R
cd ..
R CMD INSTALL SeuratRust
```

## Compare against Seurat

```r
library(Seurat)
library(SeuratRust)
library(Matrix)

mat <- Matrix::sparseMatrix(i = c(0, 2, 1), p = c(0, 1, 2, 3), x = 1:3, dims = c(3, 3))
all.equal(
  Seurat:::LogNorm(mat, 1e4, FALSE),
  SeuratRust::LogNorm(mat, 1e4, FALSE)
)
```

Parity and benchmark tests live in the parent package under
`tests/testthat/test_rust_cpp_*.R` and require `SeuratRust` in `Suggests`.

## Layout

| Path | Role |
|------|------|
| `src/rust/` | extendr crate (Rust kernels) |
| `src/cpp/` | ModularityOptimizer C++ bridge |
| `src/entrypoint.c` | Links Rust staticlib into `SeuratRust.so` |
| `R/native.R` | High-level R API matching Seurat's RcppExports |
| `R/extendr-wrappers.R` | Generated low-level `.Call` wrappers |

Seurat itself is C++/Rcpp-only; no Rust toolchain is required to build the main
package.
