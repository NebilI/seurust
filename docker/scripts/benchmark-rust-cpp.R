#!/usr/bin/env Rscript
# Print Seurat (C++) vs RSeurat timing for ported routines.
# Optional gate: SEURAT_REQUIRE_RUST_FASTER=1 fails when Rust is slower.

system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")

suppressPackageStartupMessages({
  devtools::load_all(recompile = FALSE, quiet = TRUE)
  library(RSeurat)
  library(Matrix)
})

source("tests/testthat/helper-benchmark.R", local = TRUE)

require_rust_faster <- identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")
failures <- character(0)

run_bench <- function(label, cpp_fn, rust_fn, tolerance = 0.95, ...) {
  bench <- benchmark_rust_cpp(cpp_fn = cpp_fn, rust_fn = rust_fn, ...)
  line <- format_benchmark(bench, label)
  cat(line, "\n", sep = "")
  if (require_rust_faster && bench$rust_vs_cpp < tolerance) {
    failures <<- c(failures, line)
  }
  invisible(bench)
}

run_compute_snn_bench <- function(n_cells, enforce_faster = FALSE, tolerance = 0.95, ...) {
  bench <- benchmark_compute_snn(n_cells = n_cells, ...)
  label <- attr(bench, "label")
  line <- format_benchmark(bench, label)
  cat(line, "\n", sep = "")
  if (require_rust_faster && isTRUE(enforce_faster) && bench$rust_vs_cpp < tolerance) {
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
  rust_fn = function() do.call(RSeurat::RunModularityClusteringCpp, c(list(SNN = connections), modularity_args)),
  n_warmup = 2L,
  tolerance = 0.95
)

cat("\n==> LogNorm\n")
mat <- as(matrix(1:160000, ncol = 400, nrow = 400), "sparseMatrix")
run_bench(
  "LogNorm (400x400 sparse)",
  cpp_fn = function() Seurat:::LogNorm(mat, 1e4, display_progress = FALSE),
  rust_fn = function() RSeurat::LogNorm(mat, 1e4, display_progress = FALSE)
)

cat("\n==> Dense matrix kernels\n")
set.seed(11)
dense_mat <- matrix(rnorm(2000 * 300), nrow = 2000, ncol = 300)
run_bench(
  "Standardize (2000x300 dense)",
  cpp_fn = function() Seurat:::Standardize(dense_mat, display_progress = FALSE),
  rust_fn = function() RSeurat::Standardize(dense_mat, display_progress = FALSE),
  n_warmup = 1L,
  n_reps = 25L,
  tolerance = 0.95
)
run_bench(
  "FastCov (2000x120 dense)",
  cpp_fn = function() Seurat:::FastCov(dense_mat[, 1:120], center = TRUE),
  rust_fn = function() RSeurat::FastCov(dense_mat[, 1:120], center = TRUE),
  n_warmup = 1L,
  n_reps = 25L
)
run_bench(
  "RowVar (2000x300 dense)",
  cpp_fn = function() Seurat:::RowVar(dense_mat),
  rust_fn = function() RSeurat::RowVar(dense_mat),
  n_warmup = 1L,
  n_reps = 25L,
  tolerance = 0.95
)

cat("\n==> fast_dist\n")
set.seed(12)
dist_x <- matrix(rnorm(2500 * 30), nrow = 2500, ncol = 30)
dist_y <- matrix(rnorm(2500 * 30), nrow = 2500, ncol = 30)
dist_neighbors <- replicate(2500, as.numeric(sample.int(2500, 20)), simplify = FALSE)
run_bench(
  "fast_dist (2500 cells x 30 dims, k=20)",
  cpp_fn = function() Seurat:::fast_dist(x = dist_x, y = dist_y, n = dist_neighbors),
  rust_fn = function() RSeurat::fast_dist(x = dist_x, y = dist_y, n = dist_neighbors),
  n_warmup = 1L,
  n_reps = 25L,
  tolerance = 0.95
)

cat("\n==> FastSparseRowScale\n")
scale_mat <- as(
  Matrix::rsparsematrix(nrow = 2000, ncol = 2500, density = 0.12, rand.x = stats::runif),
  "dgCMatrix"
)
run_bench(
  "FastSparseRowScale (2000x2500 sparse)",
  cpp_fn = function() {
    Seurat:::FastSparseRowScale(
      scale_mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
    )
  },
  rust_fn = function() {
    RSeurat::FastSparseRowScale(
      scale_mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
    )
  },
  n_warmup = 1L,
  tolerance = 0.95
)

cat("\n==> SparseRowVar2\n")
var_mat <- as(
  Matrix::rsparsematrix(nrow = 2000, ncol = 2500, density = 0.12, rand.x = stats::runif),
  "dgCMatrix"
)
mu <- Matrix::rowMeans(var_mat)
run_bench(
  "SparseRowVar2 (2000x2500 sparse)",
  cpp_fn = function() {
    Seurat:::SparseRowVar2(var_mat, mu = mu, display_progress = FALSE)
  },
  rust_fn = function() {
    RSeurat::SparseRowVar2(var_mat, mu = mu, display_progress = FALSE)
  },
  n_warmup = 1L,
  tolerance = 0.95
)

cat("\n==> ComputeSNN\n")
run_compute_snn_bench(500L, n_warmup = 2L, enforce_faster = FALSE)
run_compute_snn_bench(2000L, n_warmup = 1L, enforce_faster = TRUE, tolerance = 0.95)

cat("\n==> row_sum_dgcmatrix\n")
big <- sparseMatrix(
  i = sample.int(3000, 500000, replace = TRUE),
  j = sample.int(800, 500000, replace = TRUE),
  x = runif(500000),
  dims = c(3000L, 800L)
)
bx <- slot(big, "x")
bi <- slot(big, "i")
run_bench(
  "row_sum_dgcmatrix (3000x800 sparse)",
  cpp_fn = function() Seurat:::row_sum_dgcmatrix(bx, bi, nrow(big), ncol(big)),
  rust_fn = function() RSeurat::row_sum_dgcmatrix(bx, bi, nrow(big), ncol(big))
)

cat("\n==> Integration helpers\n")
set.seed(13)
n_cells <- 1200L
n_anchors <- 1800L
k_weight <- 20L
cells2 <- as.numeric(seq_len(n_cells) - 1L)
distances <- matrix(runif(n_cells * k_weight), nrow = n_cells, ncol = k_weight)
cell_index <- matrix(sample.int(n_anchors, n_cells * k_weight, replace = TRUE), nrow = n_cells)
storage.mode(cell_index) <- "double"
anchor_cells2 <- paste0("cell", sample.int(500L, n_anchors, replace = TRUE))
integration_rownames <- paste0("cell", sample.int(500L, n_anchors, replace = TRUE))
anchor_score <- runif(n_anchors, min = 0.1, max = 1)
run_bench(
  "FindWeightsC (1200 cells, 1800 anchors, min_dist=0)",
  cpp_fn = function() {
    Seurat:::FindWeightsC(
      cells2 = cells2, distances = distances, anchor_cells2 = anchor_cells2,
      integration_matrix_rownames = integration_rownames, cell_index = cell_index,
      anchor_score = anchor_score, min_dist = 0, sd = 1, display_progress = FALSE
    )
  },
  rust_fn = function() {
    RSeurat::FindWeightsC(
      cells2 = cells2, distances = distances, anchor_cells2 = anchor_cells2,
      integration_matrix_rownames = integration_rownames, cell_index = cell_index,
      anchor_score = anchor_score, min_dist = 0, sd = 1, display_progress = FALSE
    )
  },
  n_warmup = 1L,
  n_reps = 10L,
  tolerance = 0.95
)

score_nn <- matrix(sample.int(800L, 800L * 20L, replace = TRUE), nrow = 800L, ncol = 20L)
storage.mode(score_nn) <- "double"
score_snn <- Seurat:::ComputeSNN(score_nn, prune = 0.01)
query_pca <- matrix(rnorm(30L * 800L), nrow = 30L, ncol = 800L)
query_dists <- matrix(abs(rnorm(800L * 20L)), nrow = 800L, ncol = 20L)
corrected_nns <- matrix(sample.int(800L, 800L * 20L, replace = TRUE), nrow = 800L, ncol = 20L)
storage.mode(corrected_nns) <- "double"
run_bench(
  "ScoreHelper (800 cells x 30 dims, k=20)",
  cpp_fn = function() {
    Seurat:::ScoreHelper(
      snn = score_snn, query_pca = query_pca, query_dists = query_dists,
      corrected_nns = corrected_nns, k_snn = 20L, subtract_first_nn = FALSE,
      display_progress = FALSE
    )
  },
  rust_fn = function() {
    RSeurat::ScoreHelper(
      snn = score_snn, query_pca = query_pca, query_dists = query_dists,
      corrected_nns = corrected_nns, k_snn = 20L, subtract_first_nn = FALSE,
      display_progress = FALSE
    )
  },
  n_warmup = 1L,
  n_reps = 10L,
  tolerance = 0.95
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
