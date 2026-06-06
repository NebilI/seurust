#!/usr/bin/env Rscript
# Example scRNA-seq workflow using Seurat's C++/Rcpp native kernels.
#
# Mirrors a typical single-cell analysis: QC, normalization, variable features,
# scaling, PCA, SNN graph construction, Louvain clustering, UMAP, and a small
# batch-integration demo. Only the performance-critical steps call into C++.
#
# Run from repo root:
#   Rscript examples/scrna_workflow_cpp.R
#
# Or inside the dev container:
#   docker compose -f docker/docker-compose.yml run --rm rust-dev \
#     Rscript examples/scrna_workflow_cpp.R

source("examples/helpers/scrna_common.R", local = TRUE)
bootstrap_example_env()

cat("==> scRNA-seq workflow (C++ backend)\n\n")
out <- run_scrna_workflow(
  backend = make_backend("cpp"),
  output_file = "examples/output/cpp_results.rds"
)

cat("\nDone. Compare with examples/scrna_workflow_rust.R using compare_scrna_workflows.R\n")
