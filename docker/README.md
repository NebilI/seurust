# Docker development images

These images are for local seurust development. CI for pull requests uses
GitHub-hosted runners (`seurust Checks`), not these images.

| Image | Dockerfile | Purpose |
|-------|------------|---------|
| `seurust-rcpp:dev` | `Dockerfile.rcpp` | R toolchain base (compilers and common libraries) |
| `seurust-rust:dev` | `Dockerfile.rust` | Above, plus Rust stable, bindgen/clang, and `rextendr` |

## Quick start

From the repository root:

```sh
docker compose -f docker/docker-compose.yml build
docker compose -f docker/docker-compose.yml run --rm rust-dev
```

`rcpp-dev` is only the base image. `rust-dev` is the image to use for seurust.

## Install and test seurust

Inside the container, or in one shot from the host:

```sh
bash docker/scripts/build-and-test-rust.sh
```

```sh
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  bash docker/scripts/build-and-test-rust.sh
```

That installs seurust's R dependencies, runs `R CMD INSTALL seurust`, runs the
Rust library tests, and runs `seurust/tests/testthat`.

Compare kernels with an installed copy of Seurat from the repo root:

```sh
Rscript scripts/bench-all-kernels.R
```

Seurat is not built from this repository. `docker/scripts/bootstrap-dev-env.R`
installs it when a script needs it.

Refresh vendored crates if Rust dependencies changed:

```sh
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  bash docker/scripts/vendor-seurust.sh
```

## Rust crate

The crate lives in `seurust/src/rust/`. `seurust/configure` generates `Makevars`
and runs `cargo build`. Development installs also run the `document`
binary to refresh `R/extendr-wrappers.R`.

```sh
cargo test --manifest-path seurust/src/rust/Cargo.toml --lib
```

```sh
cd seurust/src/rust && cargo run --bin document --release && cd ../../..
```

## Build without Compose

```sh
docker build -f docker/Dockerfile.rcpp -t seurust-rcpp:dev .
docker build -f docker/Dockerfile.rust -t seurust-rust:dev .
docker run --rm -it -v "$(pwd):/workspace" -w /workspace seurust-rust:dev
```

On Windows PowerShell, replace `$(pwd)` with `${PWD}`.

## Notes

- The base image is `rocker/r2u:jammy`.
- The rust image mounts a named volume at `seurust/src/rust/target` so Cargo artifacts stay off the bind mount.
- `seurust/src/rust/vendor.tar.xz` is the offline crate source used by source builds. Do not delete it.
