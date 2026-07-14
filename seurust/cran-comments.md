## Test environments

* local via Docker: `docker compose -f docker/docker-compose.yml run --rm seurust-cran`
  (Ubuntu 22.04, R 4.6.1, rustc stable)
* GitHub Actions (ubuntu-latest), R release — `build-seurust-cran.yaml`

## R CMD check results

There were no ERRORs.

Notes / warnings expected for this package:

* Compiled code uses Rust (extendr) plus a small C++ bridge for
  ModularityOptimizer. `SystemRequirements` lists Cargo and rustc (>= 1.81).
* Source tarball includes vendored Rust crates in
  `src/rust/vendor.tar.xz` (~2 MB) so CRAN can build offline.
* `src/Makevars` / `src/Makevars.in` use `.NOTPARALLEL` (GNU make) to avoid
  races between cargo `build.rs` and cleanup. This matches the extendr
  CRAN template and is required for reliable parallel `make`.

## Downstream dependencies

There are currently no reverse dependencies on CRAN.
