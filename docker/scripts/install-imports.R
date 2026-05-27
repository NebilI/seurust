#!/usr/bin/env Rscript
# Install only Depends/Imports/LinkingTo (enough for compile_dll on a dev checkout).
repos <- c(
  CRAN = "https://cloud.r-project.org",
  satijalab = "https://satijalab.r-universe.dev"
)
options(repos = repos)
pkg_root <- Sys.getenv("SEURAT_PKG_ROOT", unset = "/workspace")
remotes::install_deps(
  pkg_root,
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade = "never"
)
message("Imports installed.")
