# seurust

Rust/extendr backend for Seurat's performance-critical native routines. Install
alongside [Seurat](https://satijalab.org/seurat) to use the same function
signatures with a Rust backend.

## Requirements

- R (>= 4.0.0)
- Rust toolchain: [rustc](https://rust-lang.org/tools/install/) and Cargo (>= 1.81)
- On Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) plus Rust

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
install.packages("path/to/seurust_0.1.2.tar.gz", repos = NULL, type = "source")
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

## Further documentation

Development, packaging, and contribution notes live in the
[GitHub repository](https://github.com/NebilI/seurust).
