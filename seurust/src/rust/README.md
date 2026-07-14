# seurust (Rust crate)

Rust native kernels for the [seurust](https://github.com/NebilI/seurust/tree/main/seurust) R package.

## Install the R package

From [CRAN](https://cran.r-project.org/package=seurust):

```r
install.packages("seurust")
```

Or from r-universe:

```r
install.packages(
  "seurust",
  repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
)
```

## crates.io

This crate is published to [crates.io](https://crates.io/crates/seurust) on GitHub
Releases via `.github/workflows/publish-seurust-crate.yaml`. It is primarily
consumed as an R extension (`crate-type = ["rlib", "staticlib"]`), not as a
standalone Rust library API.
