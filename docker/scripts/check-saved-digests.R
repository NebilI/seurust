#!/usr/bin/env Rscript
system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")
suppressPackageStartupMessages({
  devtools::load_all(recompile = FALSE, quiet = TRUE)
  library(RSeurat)
})
cpp <- readRDS("examples/output/cpp_results.rds")
rust <- readRDS("examples/output/rust_results.rds")
cat("cluster digest match:", identical(cpp$cluster_digest, rust$cluster_digest), "\n")
cat("snn digest match:", identical(cpp$snn_digest, rust$snn_digest), "\n")
cat("umap digest match:", identical(cpp$umap_digest, rust$umap_digest), "\n")
