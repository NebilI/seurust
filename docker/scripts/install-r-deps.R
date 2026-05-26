#!/usr/bin/env Rscript
# Install R packages needed to compile Seurat's Rcpp/C++ sources from a dev checkout.
# Full package checks need additional Imports; run install-r-deps.R --full for those.

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
  "rcmdcheck",
  "withr",
  "cli"
)
install.packages(build_pkgs, Ncpus = max(1L, parallel::detectCores() - 1L))

message("Installing SeuratObject from satijalab r-universe...")
install.packages("SeuratObject", repos = "https://satijalab.r-universe.dev")

if (full) {
  message("Installing all package dependencies from DESCRIPTION (this can take a while)...")
  pkg_root <- Sys.getenv("SEURAT_PKG_ROOT", unset = "/workspace")
  if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    stop("SEURAT_PKG_ROOT does not contain DESCRIPTION: ", pkg_root)
  }
  remotes::install_deps(pkg_root, dependencies = TRUE, upgrade = "never")
}

message("Done.")
