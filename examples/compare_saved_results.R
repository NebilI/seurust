#!/usr/bin/env Rscript
# Compare timing and fingerprints from examples/output/{cpp,rust}_results.rds
# without re-running the workflows.

cpp_file <- "examples/output/cpp_results.rds"
rust_file <- "examples/output/rust_results.rds"
if (!file.exists(cpp_file) || !file.exists(rust_file)) {
  stop("Run scrna_workflow_cpp.R and scrna_workflow_rust.R first.")
}

cpp <- readRDS(cpp_file)
rust <- readRDS(rust_file)

compare_fields <- c(
  "n_cells", "n_clusters", "cluster_table", "cluster_digest",
  "norm_digest", "snn_digest", "snn_nnz", "umap_digest"
)

cat("==> Output parity\n")
all_ok <- TRUE
for (field in compare_fields) {
  ok <- identical(cpp[[field]], rust[[field]])
  all_ok <- all_ok && ok
  cat(sprintf("  %-20s %s\n", field, if (ok) "OK" else "MISMATCH"))
}

integration_ok <- identical(cpp$integration, rust$integration)
all_ok <- all_ok && integration_ok
cat(sprintf("  %-20s %s\n", "integration", if (integration_ok) "OK" else "MISMATCH"))

cat("\n==> Timing comparison (seconds)\n")
steps <- names(cpp$timings)
cat(sprintf("%-28s %10s %10s %10s\n", "Step", "C++", "Rust", "Rust/C++"))
cat(strrep("-", 62), "\n", sep = "")
for (step in steps) {
  t_cpp <- cpp$timings[[step]]
  t_rust <- rust$timings[[step]]
  cat(sprintf(
    "%-28s %10.3f %10.3f %10.2fx\n",
    step, t_cpp, t_rust, t_rust / t_cpp
  ))
}
cat(strrep("-", 62), "\n", sep = "")
cat(sprintf(
  "%-28s %10.3f %10.3f %10.2fx\n",
  "Total",
  cpp$total_native_seconds,
  rust$total_native_seconds,
  rust$total_native_seconds / cpp$total_native_seconds
))

if (!all_ok) {
  stop("Output mismatch between C++ and Rust workflows.", call. = FALSE)
}
cat("\nAll outputs match.\n")
