## seurust 0.1.0 (version line restart)

This submission replaces the earlier 0.1.x release candidates. Prior GitHub
release tags `v0.1.2`, `v0.1.4`, and `v0.1.5` were removed; the canonical
release is **`v0.1.0`** on GitHub and r-universe.

## Test environments

* GitHub Actions (ubuntu-latest), R release — `build-seurust-cran.yaml`, merge checks
* Local Docker CRAN check via `docker compose … run --rm seurust-cran`

## R CMD check results

There were no ERRORs.

Notes expected for this package:

* New submission / compiled code
* SystemRequirements for Rust toolchain on source installs
