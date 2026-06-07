#!/usr/bin/env Rscript
system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")
suppressPackageStartupMessages({
  devtools::load_all(recompile = FALSE, quiet = TRUE)
  library(RSeurat)
  library(Matrix)
})
source("examples/helpers/scrna_common.R", local = TRUE)
set.seed(42)
counts <- simulate_scrna_counts(2500, 2000, 42)
cpp_norm <- Seurat:::LogNorm(counts, 1e4, FALSE)
rust_norm <- RSeurat::LogNorm(counts, 1e4, FALSE)
cat("lognorm max diff:", max(abs(cpp_norm - rust_norm)), "\n")
cat("lognorm all.equal:", isTRUE(all.equal(as.matrix(cpp_norm), as.matrix(rust_norm))), "\n")
