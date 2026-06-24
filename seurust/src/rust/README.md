# seurust

Rust native kernels for single-cell analysis, powering the [seurust](https://github.com/NebilI/seurust/tree/main/seurust) R package.

This crate implements performance-critical routines used by [Seurat](https://satijalab.org/seurat): sparse matrix normalization and scaling, shared-nearest-neighbor graph construction, modularity clustering, and batch-integration helpers.

## Usage

The crate is primarily consumed as an R extension via [extendr](https://github.com/extendr/extendr). Install the R package for the supported workflow:

```r
remotes::install_github("NebilI/seurust", subdir = "seurust")
```

Building from source requires R, RcppEigen, and a Rust toolchain (see the seurust package documentation).

## License

MIT
