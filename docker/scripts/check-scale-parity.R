#!/usr/bin/env Rscript
system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")
source("examples/helpers/scrna_common.R", local = TRUE)
bootstrap_example_env()
set.seed(42)
counts <- simulate_scrna_counts(2500, 2000, 42)
cpp_norm <- Seurat:::LogNorm(counts, 1e4, FALSE)
rust_norm <- seurust::LogNorm(counts, 1e4, FALSE)
hvf <- rownames(cpp_norm)[1:2000]
cpp_scale_on_cpp <- Seurat:::FastSparseRowScale(cpp_norm[hvf, , drop = FALSE], TRUE, TRUE, 10, FALSE)
rust_scale_on_cpp <- seurust::FastSparseRowScale(cpp_norm[hvf, , drop = FALSE], TRUE, TRUE, 10, FALSE)
rust_scale_on_rust <- seurust::FastSparseRowScale(rust_norm[hvf, , drop = FALSE], TRUE, TRUE, 10, FALSE)
cat("rust scale on cpp norm max diff:", max(abs(cpp_scale_on_cpp - rust_scale_on_cpp)), "\n")
cat("rust vs cpp scale on own norms max diff:", max(abs(cpp_scale_on_cpp - rust_scale_on_rust)), "\n")
