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
install.packages("path/to/seurust_0.2.0.tar.gz", repos = NULL, type = "source")
```

## Example

```r
library(Seurat)
library(seurust)
library(Matrix)

mat <- Matrix::sparseMatrix(i = c(1, 3, 2), j = c(1, 2, 3), x = 1:3, dims = c(3, 3))
# Seurat 5.6.0's C++ LogNorm fills the x slot; seurust still returns a dgCMatrix.
mat_seurat <- mat
mat_seurat@x <- Seurat:::LogNorm(mat@x, mat@p, 1e4, 1L, FALSE)
all.equal(mat_seurat, seurust::LogNorm(mat, 1e4, FALSE))
```

## Differences from Seurat's C++ entry points

Most exported functions take the same arguments and return the same kind of
object as the Seurat 5.6.0 function of the same name. Integer inputs are accepted
wherever Seurat's Rcpp signature would coerce them (for example integer
neighbor ids in `fast_dist` or an integer `nn_ranked` matrix in `ComputeSNN`).

Deliberate differences:

- Seurat 5.6.0's C++ `LogNorm` takes CSC `x`/`p` slots and returns the
  normalized non-zeros. `seurust::LogNorm` still accepts a `dgCMatrix` and
  returns a `dgCMatrix` (the same object `Seurat::LogNormalize` builds).
- Seurat's `RowMergeMatrices` only accepts row-compressed (`RsparseMatrix`)
  inputs and errors on a `dgCMatrix`. `seurust::RowMergeMatrices` accepts
  either layout and converts first; for `RsparseMatrix` inputs the result is
  identical.
- `AnnoySearchCpp` is not ported.

## Further documentation

Development, packaging, and contribution notes live in the
[GitHub repository](https://github.com/NebilI/seurust).
