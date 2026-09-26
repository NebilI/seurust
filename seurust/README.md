# seurust

Rust/extendr backend for Seurat's performance-critical native routines. Install
alongside [Seurat](https://satijalab.org/seurat) to use the same function
signatures with a Rust backend.

## Requirements

- R (>= 4.0.0)
- Rust toolchain: [rustc](https://rust-lang.org/tools/install/) and Cargo (>= 1.81)
- On Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) plus Rust

No C++ compiler or `RcppEigen` headers are needed; all compiled code is Rust.

## Install

### From CRAN (after acceptance)

```r
install.packages("seurust")
```

### From r-universe

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

Download `seurust_*.tar.gz` from
[GitHub Releases](https://github.com/NebilI/seurust/releases), then:

```r
install.packages("path/to/seurust_0.1.3.tar.gz", repos = NULL, type = "source")
```

## Example

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

## Differences from Seurat's C++ entry points

Every exported function takes the same arguments and returns the same kind of
object as the Seurat function of the same name. Integer inputs are accepted
wherever Seurat's Rcpp signature would coerce them (for example integer
neighbor ids in `fast_dist` or an integer `nn_ranked` matrix in `ComputeSNN`).

One deliberate difference: Seurat's `RowMergeMatrices` only accepts
row-compressed (`RsparseMatrix`) inputs and errors on a `dgCMatrix`.
`seurust::RowMergeMatrices` accepts either layout and converts first; for
`RsparseMatrix` inputs the result is identical.

## Further documentation

Development, packaging, and contribution notes live in the
[GitHub repository](https://github.com/NebilI/seurust).
