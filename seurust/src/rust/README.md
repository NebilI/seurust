# seurust (Rust crate)

Rust native kernels for the [seurust](https://github.com/NebilI/seurust/tree/main/seurust) R package.

## Install the R package

From [r-universe](https://NebilI.r-universe.dev/seurust):

```r
install.packages(
  "seurust",
  repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
)
```

## crates.io

This crate is published to [crates.io](https://crates.io/crates/seurust) on GitHub
Releases via `.github/workflows/release-seurust.yaml`. It is primarily
consumed as an R extension (`crate-type = ["rlib", "staticlib"]`), not as a
standalone Rust library API.
