#!/usr/bin/env Rscript
# Print Seurat (C++) vs SeuratRust timing for ported routines.
# Optional gate: SEURAT_REQUIRE_RUST_FASTER=1 fails when Rust is slower.

system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")

suppressPackageStartupMessages({
  devtools::load_all(recompile = FALSE, quiet = TRUE)
  library(SeuratRust)
  library(Matrix)
})

source("tests/testthat/helper-benchmark.R", local = TRUE)

require_rust_faster <- identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")
failures <- character(0)

run_bench <- function(label, cpp_fn, rust_fn, ...) {
  bench <- benchmark_rust_cpp(cpp_fn = cpp_fn, rust_fn = rust_fn, ...)
  line <- format_benchmark(bench, label)
  cat(line, "\n", sep = "")
  if (require_rust_faster && bench$rust_vs_cpp < 1) {
    failures <<- c(failures, line)
  }
  invisible(bench)
}

cat("==> Modularity clustering\n")
node1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
           1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 8, 8, 8, 9, 13,
           14, 14, 15, 15, 18, 18, 19, 20, 20, 22, 22, 23, 23, 23, 23, 23, 24,
           24, 24, 25, 26, 26, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32)
node2 <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 17, 19, 21, 31, 2, 3, 7, 13,
           17, 19, 21, 30, 3, 7, 8, 9, 13, 27, 28, 32, 7, 12, 13, 6, 10, 6, 10,
           16, 16, 30, 32, 33, 33, 33, 32, 33, 32, 33, 32, 33, 33, 32, 33, 32,
           33, 25, 27, 29, 32, 33, 25, 27, 31, 31, 29, 33, 33, 31, 33, 32, 33,
           32, 33, 32, 33, 33)
connections <- sparseMatrix(i = node2 + 1, j = node1 + 1, x = 1.0)
modularity_args <- list(
  modularityFunction = 1L,
  resolution = 1.0,
  algorithm = 3L,
  nRandomStarts = 5L,
  nIterations = 50L,
  randomSeed = 42L,
  printOutput = FALSE,
  edgefilename = ""
)
run_bench(
  "Modularity (alg 3, 5 starts x 50 iters)",
  cpp_fn = function() do.call(Seurat:::RunModularityClusteringCpp, c(list(SNN = connections), modularity_args)),
  rust_fn = function() do.call(SeuratRust::RunModularityClusteringCpp, c(list(SNN = connections), modularity_args)),
  n_warmup = 2L,
  n_reps = 10L
)

cat("\n==> LogNorm\n")
mat <- as(matrix(1:160000, ncol = 400, nrow = 400), "sparseMatrix")
run_bench(
  "LogNorm (400x400 sparse)",
  cpp_fn = function() Seurat:::LogNorm(mat, 1e4, display_progress = FALSE),
  rust_fn = function() SeuratRust::LogNorm(mat, 1e4, display_progress = FALSE)
)

cat("\n==> ComputeSNN\n")
set.seed(1)
nn <- matrix(
  sample.int(500, 500 * 20, replace = TRUE),
  nrow = 500,
  ncol = 20
)
storage.mode(nn) <- "double"
run_bench(
  "ComputeSNN (500 cells, k=20)",
  cpp_fn = function() Seurat:::ComputeSNN(nn, 0.01),
  rust_fn = function() SeuratRust::ComputeSNN(nn, 0.01)
)

cat("\n==> row_sum_dgcmatrix\n")
big <- sparseMatrix(
  i = sample.int(3000, 50000, replace = TRUE),
  j = sample.int(800, 50000, replace = TRUE),
  x = runif(50000),
  dims = c(3000L, 800L)
)
bx <- slot(big, "x")
bi <- slot(big, "i")
run_bench(
  "row_sum_dgcmatrix (3000x800 sparse)",
  cpp_fn = function() Seurat:::row_sum_dgcmatrix(bx, bi, nrow(big), ncol(big)),
  rust_fn = function() SeuratRust::row_sum_dgcmatrix(bx, bi, nrow(big), ncol(big))
)

if (length(failures) > 0) {
  stop(
    "Rust was slower than C++ for:\n",
    paste0("  - ", failures, collapse = "\n"),
    "\nSet SEURAT_REQUIRE_RUST_FASTER=0 to report without failing."
  )
}

cat("\nBenchmark complete.\n")
cat("Ratio > 1.0 means Rust is faster. Set SEURAT_REQUIRE_RUST_FASTER=1 to enforce.\n")
