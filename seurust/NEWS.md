# seurust 0.1.2

* Update maintainer contact for CRAN submissions.
* Link the Rust static library with Apple ld (`-force_load`) so macOS
  source builds succeed.
* Carry forward the 0.1.1 packaging fixes for Windows configure and portable
  Makevars.

# seurust 0.1.1

* Fix Windows source installs for CRAN by adding `configure.win` /
  `cleanup.win` so `Makevars.win` is generated and the Rust static library
  is linked.
* Remove non-portable `.NOTPARALLEL` from Makevars templates.
* Clarify DESCRIPTION / README wording for CRAN incoming checks.

# seurust 0.1.0

* Initial public release of Rust/extendr native kernels matching Seurat's
  internal C++ API (`LogNorm`, `FastSparseRowScale`, `ComputeSNN`,
  `IntegrateDataC`, `RunModularityClusteringCpp`, and related helpers).
* Offline CRAN builds via vendored Rust crates (`src/rust/vendor.tar.xz`).
* Installable from GitHub Releases, r-universe, and (after acceptance) CRAN.
