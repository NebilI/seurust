#!/usr/bin/env Rscript
# Run both scRNA-seq workflow examples and verify identical outputs + timing.
#
# Run from repo root:
#   Rscript examples/compare_scrna_workflows.R
#
# Or inside the dev container:
#   docker compose -f docker/docker-compose.yml run --rm rust-dev \
#     Rscript examples/compare_scrna_workflows.R

source("examples/helpers/scrna_common.R", local = TRUE)
bootstrap_example_env()

out_dir <- "examples/output"
cpp_file <- file.path(out_dir, "cpp_results.rds")
rust_file <- file.path(out_dir, "rust_results.rds")

run_script <- function(path) {
  status <- system2("Rscript", path, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Script failed: ", path, call. = FALSE)
  }
}

cat("==> Running C++ workflow...\n\n")
run_script("examples/scrna_workflow_cpp.R")

cat("\n==> Running Rust workflow...\n\n")
run_script("examples/scrna_workflow_rust.R")

cpp <- readRDS(cpp_file)
rust <- readRDS(rust_file)

compare_fields <- c(
  "n_cells", "n_clusters", "cluster_table", "cluster_digest",
  "norm_digest", "snn_digest", "snn_nnz", "umap_digest"
)

cat("\n==> Output parity\n")
all_ok <- TRUE
for (field in compare_fields) {
  ok <- identical(cpp[[field]], rust[[field]])
  all_ok <- all_ok && ok
  status <- if (ok) "OK" else "MISMATCH"
  cat(sprintf("  %-20s %s\n", field, status))
  if (!ok) {
    cat("    C++ :", cpp[[field]], "\n")
    cat("    Rust:", rust[[field]], "\n")
  }
}

integration_ok <- identical(
  cpp$integration$integrated_digest,
  rust$integration$integrated_digest
) && identical(
  cpp$integration$weights_digest,
  rust$integration$weights_digest
)
all_ok <- all_ok && integration_ok
cat(sprintf("  %-20s %s\n", "integration_digest",
            if (integration_ok) "OK" else "MISMATCH"))

cat("\n==> Timing comparison (seconds)\n")
steps <- names(cpp$timings)
cat(sprintf("%-28s %10s %10s %10s\n", "Step", "C++", "Rust", "Rust/C++"))
cat(strrep("-", 62), "\n", sep = "")
for (step in steps) {
  t_cpp <- cpp$timings[[step]]
  t_rust <- rust$timings[[step]]
  ratio <- if (t_cpp > 0) t_rust / t_cpp else NA_real_
  cat(sprintf(
    "%-28s %10.3f %10.3f %10.2fx\n",
    step, t_cpp, t_rust, ratio
  ))
}
cat(strrep("-", 62), "\n", sep = "")
total_cpp <- cpp$total_native_seconds
total_rust <- rust$total_native_seconds
cat(sprintf(
  "%-28s %10.3f %10.3f %10.2fx\n",
  "Total",
  total_cpp,
  total_rust,
  total_rust / total_cpp
))

if (!all_ok) {
  stop("Output mismatch between C++ and Rust workflows.", call. = FALSE)
}

cat("\nAll outputs match. Ratio < 1.0 means Rust was faster for that step.\n")
