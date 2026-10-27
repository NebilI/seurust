#!/usr/bin/env Rscript
# Install Depends/Imports/LinkingTo for the seurust package in this repository.
repos <- c(
  CRAN = "https://cloud.r-project.org",
  satijalab = "https://satijalab.r-universe.dev"
)
options(repos = repos)
pkg <- "/workspace/seurust"
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_deps(
  pkg,
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade = "never"
)
message("seurust Imports installed.")
