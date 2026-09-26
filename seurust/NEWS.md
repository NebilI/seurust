# seurust 0.1.3

## Seurat compatibility fixes

* `fast_dist()` accepts integer neighbor ids (previously it panicked unless
  ids were doubles), keeps the names of the neighbor list like Seurat, and
  reports out-of-range or `NA` ids as R errors.
* `GraphToNeighborHelper()` now reads each cell's neighbors from its row, as
  Seurat does; it previously read columns, which gave different neighbor sets
  for asymmetric graphs. Ragged rows raise an R error instead of a panic.
* `DirectSNNToFile()` returns the SNN graph as a `dgCMatrix` (it returned the
  raw `x`/`i`/`p`/`Dim` list) and validates `nn_ranked` like `ComputeSNN()`.
* `DirectSNNToFile()` and `WriteEdgeFile()` format edge weights like C++
  `setprecision(15)`, so edge files are byte-identical to Seurat's. File
  errors are R errors instead of panics.
* Integer matrices and vectors are coerced to double before reaching Rust
  wherever Seurat's Rcpp signature would coerce them (`ComputeSNN()`,
  `FastCov()`, `SparseRowVar2()` means, and others).
* `RowMergeMatrices()` still accepts a `dgCMatrix` (Seurat requires an
  `RsparseMatrix`); this is now documented.

## Performance

* `FastCov()` and `FastCovMats()` read R memory directly and use a blocked,
  multi-threaded matrix product (`FastCov()` only computes the upper
  triangle). They were several times slower than Seurat and are now faster.
* `FastRBind()` copies columns straight into the result instead of going
  through element-wise `ndarray` conversions (about 50x faster than 0.1.2
  on large inputs).

## Build

* Removed the unused Eigen C++ bridge. Source builds no longer need a C++
  compiler or `RcppEigen`, which fixes builds where `c++` is Clang without a
  usable C++ standard library (`'cmath' file not found`).

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
