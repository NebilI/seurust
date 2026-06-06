#!/usr/bin/env Rscript
# Example scRNA-seq workflow using RSeurat's Rust/extendr native kernels.
#
# Identical analysis steps to scrna_workflow_cpp.R; only the backend namespace
# differs (RSeurat:: vs Seurat:::). Outputs should match the C++ run bit-for-bit
# on ported routines.
#
# Run from repo root:
#   Rscript examples/scrna_workflow_rust.R
#
# Or inside the dev container:
#   docker compose -f docker/docker-compose.yml run --rm rust-dev \
#     Rscript examples/scrna_workflow_rust.R

source("examples/helpers/scrna_common.R", local = TRUE)
bootstrap_example_env()

cat("==> scRNA-seq workflow (Rust backend)\n\n")
out <- run_scrna_workflow(
  backend = make_backend("rust"),
  output_file = "examples/output/rust_results.rds"
)

cat("\nDone. Compare with examples/scrna_workflow_cpp.R using compare_scrna_workflows.R\n")
