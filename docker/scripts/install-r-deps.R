#!/usr/bin/env Rscript
# Install R packages used by the seurust development images.
# Pass --full to also install Seurat for parity checks.

args <- commandArgs(trailingOnly = TRUE)
full <- "--full" %in% args

repos <- c(
  CRAN = "https://cloud.r-project.org",
  satijalab = "https://satijalab.r-universe.dev"
)
options(repos = repos)

message("Installing R build tooling and LinkingTo dependencies...")
build_pkgs <- c(
  "Rcpp",
  "RcppEigen",
  "RcppProgress",
  "Matrix",
  "methods",
  "remotes",
  "pkgbuild",
  "devtools",
  "roxygen2",
  "testthat",
  "microbenchmark",
  "rcmdcheck",
  "withr",
  "cli"
)
install.packages(build_pkgs, Ncpus = max(1L, parallel::detectCores() - 1L))

message("Installing SeuratObject from satijalab r-universe...")
install.packages("SeuratObject", repos = "https://satijalab.r-universe.dev")

if (full) {
  message("Installing Seurat (this repository does not vendor Seurat sources)...")
  install.packages("Seurat", Ncpus = max(1L, parallel::detectCores() - 1L))
}

message("Done.")
