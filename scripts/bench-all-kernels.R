#!/usr/bin/env Rscript
# Time every ported seurust kernel against Seurat's C++ implementation.
# Prints a markdown table: Seurat C++, seurust Rust, and speedup (C++ / Rust).

.libPaths(c(
  Sys.getenv("SEURUST_LIB", unset = "/tmp/lib-seurust"),
  Sys.getenv("R_LIBS_USER", unset = path.expand("~/R/library")),
  .libPaths()
))

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(seurust)
  library(microbenchmark)
})

seurat_fn <- function(name) getFromNamespace(name, "Seurat")

counts_sparse <- function(nrow, ncol, density = 0.05, seed = 1L) {
  set.seed(seed)
  as(
    rsparsematrix(nrow, ncol, density = density, rand.x = function(n) rpois(n, 4) + 1),
    "dgCMatrix"
  )
}

dense_mat <- function(nrow, ncol, seed = 1L) {
  set.seed(seed)
  matrix(rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
}

regular_graph <- function(n, k) {
  offs <- c(-seq_len(k), seq_len(k))
  rows <- rep(seq_len(n), each = length(offs))
  cols <- ((rows - 1L + rep(offs, times = n)) %% n) + 1L
  vals <- rep(abs(as.numeric(offs)), times = n)
  sparseMatrix(i = rows, j = cols, x = vals, dims = c(n, n))
}

nn_index <- function(n, k, seed = 1L) {
  set.seed(seed)
  nn <- matrix(sample.int(n, n * k, replace = TRUE), nrow = n, ncol = k)
  storage.mode(nn) <- "double"
  nn
}

fmt_ms <- function(us) {
  ms <- us / 1000
  if (ms >= 100) sprintf("%.0f ms", ms) else if (ms >= 10) sprintf("%.1f ms", ms) else sprintf("%.2f ms", ms)
}

time_pair <- function(cpp_fn, rust_fn, times = 7L) {
  force(cpp_fn())
  force(rust_fn())
  gc(verbose = FALSE)
  cpp <- microbenchmark(cpp_fn(), times = times, warmup = 0L)
  rust <- microbenchmark(rust_fn(), times = times, warmup = 0L)
  cpp_us <- median(cpp$time) / 1000
  rust_us <- median(rust$time) / 1000
  speedup <- cpp_us / rust_us
  list(cpp_us = cpp_us, rust_us = rust_us, speedup = speedup, times = times)
}

cat("Seurat ", as.character(packageVersion("Seurat")),
    "  seurust ", as.character(packageVersion("seurust")),
    "  R ", R.version.string, "\n\n", sep = "")

# Shared inputs. Built once so setup is not part of the timed call.
set.seed(1)
sparse <- counts_sparse(5000L, 2000L, density = 0.05, seed = 1L)
mu <- as.numeric(Matrix::rowMeans(sparse))
sigma <- sqrt(pmax(mu, 1e-8))
dense <- dense_mat(2000L, 120L, seed = 2L)
cov_a <- dense[, seq_len(80L), drop = FALSE]
cov_b <- dense_mat(2000L, 40L, seed = 3L)
rbind_a <- dense_mat(4000L, 150L, seed = 4L)
graph <- regular_graph(2000L, 15L)
nn <- nn_index(1500L, 21L, seed = 5L)
nn_file <- nn_index(400L, 16L, seed = 6L)
snn_file <- seurat_fn("ComputeSNN")(nn_file, 1 / 15)
umi <- counts_sparse(2000L, 1500L, density = 0.08, seed = 7L)
umi_per_cell <- as.numeric(rep(40, ncol(umi)))

merge_a <- as(counts_sparse(1500L, 400L, 0.04, seed = 8L), "RsparseMatrix")
merge_b <- as(counts_sparse(1200L, 350L, 0.04, seed = 9L), "RsparseMatrix")
rn1 <- paste0("g", seq_len(nrow(merge_a)))
rn2 <- paste0("g", seq(200L, length.out = nrow(merge_b)))
all_rn <- union(rn1, rn2)

replace_mat <- counts_sparse(2000L, 600L, 0.04, seed = 10L)
replace_src <- counts_sparse(2000L, 20L, 0.08, seed = 11L)
replace_idx <- as.numeric(seq(0L, 19L) * 20L)

expr <- counts_sparse(800L, 200L, 0.08, seed = 12L)
im <- counts_sparse(60L, 200L, 0.2, seed = 13L)
weights <- counts_sparse(60L, 800L, 0.05, seed = 14L)

n_cells_w <- 400L
n_anchors <- 600L
k_weight <- 10L
cells2 <- as.numeric(seq_len(n_cells_w) - 1L)
distances <- matrix(runif(n_cells_w * k_weight), nrow = n_cells_w)
cell_index <- matrix(sample.int(n_anchors, n_cells_w * k_weight, replace = TRUE), nrow = n_cells_w)
storage.mode(cell_index) <- "double"
anchor_cells2 <- paste0("cell", sample.int(80L, n_anchors, replace = TRUE))
integration_rownames <- paste0("cell", sample.int(80L, n_anchors, replace = TRUE))
anchor_score <- runif(n_anchors, min = 0.1, max = 1)

score_nn <- nn_index(300L, 12L, seed = 15L)
score_snn <- seurat_fn("ComputeSNN")(score_nn, 0)
query_pca <- dense_mat(16L, 300L, seed = 16L)
query_dists <- abs(dense_mat(300L, 12L, seed = 17L))
corrected_nns <- nn_index(300L, 12L, seed = 18L)

dist_x <- dense_mat(800L, 20L, seed = 19L)
dist_y <- dense_mat(800L, 20L, seed = 20L)
dist_n <- replicate(800L, as.numeric(sample.int(800L, 10L)), simplify = FALSE)

row_mat <- counts_sparse(6000L, 1500L, 0.04, seed = 21L)
row_x <- row_mat@x
row_i <- row_mat@i
row_nr <- nrow(row_mat)
row_nc <- ncol(row_mat)

mod_nn <- nn_index(250L, 15L, seed = 22L)
mod_snn <- seurat_fn("ComputeSNN")(mod_nn, 0.01)
embed <- dense_mat(nrow(snn_file), 10L, seed = 23L)
nearest_dist <- abs(rnorm(nrow(snn_file)))
nearest_dist[c(1L, 5L)] <- 0

edge_cpp <- tempfile(fileext = ".txt")
edge_rust <- tempfile(fileext = ".txt")
direct_cpp <- tempfile(fileext = ".txt")
direct_rust <- tempfile(fileext = ".txt")
on.exit(unlink(c(edge_cpp, edge_rust, direct_cpp, direct_rust)), add = TRUE)

cases <- list(
  list("LogNorm", "5,000 × 2,000 sparse, 5% nonzero", 5L,
    function() seurat_fn("LogNorm")(sparse, 1e4, FALSE),
    function() seurust::LogNorm(sparse, 1e4, FALSE)),
  list("Standardize", "2,000 × 120 dense", 7L,
    function() seurat_fn("Standardize")(dense, FALSE),
    function() seurust::Standardize(dense, FALSE)),
  list("FastCov", "2,000 × 80 dense", 5L,
    function() seurat_fn("FastCov")(cov_a, TRUE),
    function() seurust::FastCov(cov_a, TRUE)),
  list("FastCovMats", "2,000 × 80 and 2,000 × 40 dense", 5L,
    function() seurat_fn("FastCovMats")(cov_a, cov_b, TRUE),
    function() seurust::FastCovMats(cov_a, cov_b, TRUE)),
  list("FastRBind", "two 4,000 × 150 dense", 7L,
    function() seurat_fn("FastRBind")(rbind_a, rbind_a),
    function() seurust::FastRBind(rbind_a, rbind_a)),
  list("RowVar", "2,000 × 120 dense", 7L,
    function() seurat_fn("RowVar")(dense),
    function() seurust::RowVar(dense)),
  list("FastExpMean", "5,000 × 2,000 sparse, 5% nonzero", 5L,
    function() seurat_fn("FastExpMean")(sparse, FALSE),
    function() seurust::FastExpMean(sparse, FALSE)),
  list("SparseRowVar", "5,000 × 2,000 sparse, 5% nonzero", 5L,
    function() seurat_fn("SparseRowVar")(sparse, FALSE),
    function() seurust::SparseRowVar(sparse, FALSE)),
  list("SparseRowVar2", "5,000 × 2,000 sparse, 5% nonzero", 5L,
    function() seurat_fn("SparseRowVar2")(sparse, mu, FALSE),
    function() seurust::SparseRowVar2(sparse, mu, FALSE)),
  list("SparseRowVarStd", "5,000 × 2,000 sparse, 5% nonzero", 5L,
    function() seurat_fn("SparseRowVarStd")(sparse, mu, sigma, 10, FALSE),
    function() seurust::SparseRowVarStd(sparse, mu, sigma, 10, FALSE)),
  list("FastLogVMR", "5,000 × 2,000 sparse, 5% nonzero", 5L,
    function() seurat_fn("FastLogVMR")(sparse, FALSE),
    function() seurust::FastLogVMR(sparse, FALSE)),
  list("FastSparseRowScale", "5,000 × 2,000 sparse, 5% nonzero", 3L,
    function() seurat_fn("FastSparseRowScale")(sparse, TRUE, TRUE, 10, FALSE),
    function() seurust::FastSparseRowScale(sparse, TRUE, TRUE, 10, FALSE)),
  list("FastSparseRowScaleWithKnownStats", "5,000 × 2,000 sparse, 5% nonzero", 3L,
    function() seurat_fn("FastSparseRowScaleWithKnownStats")(sparse, mu, sigma, TRUE, TRUE, 10, FALSE),
    function() seurust::FastSparseRowScaleWithKnownStats(sparse, mu, sigma, TRUE, TRUE, 10, FALSE)),
  list("RowMergeMatrices", "1,500 × 400 and 1,200 × 350 CSR", 5L,
    function() seurat_fn("RowMergeMatrices")(merge_a, merge_b, rn1, rn2, all_rn),
    function() seurust::RowMergeMatrices(merge_a, merge_b, rn1, rn2, all_rn)),
  list("ReplaceColsC", "2,000 × 600 sparse, replace 20 columns", 7L,
    function() seurat_fn("ReplaceColsC")(replace_mat, replace_idx, replace_src),
    function() seurust::ReplaceColsC(replace_mat, replace_idx, replace_src)),
  list("GraphToNeighborHelper", "2,000 cells, 30 neighbors, symmetric", 7L,
    function() seurat_fn("GraphToNeighborHelper")(graph),
    function() seurust::GraphToNeighborHelper(graph)),
  list("RunUMISampling", "2,000 × 1,500 counts, 8% nonzero", 5L,
    function() seurat_fn("RunUMISampling")(umi, 50, FALSE, FALSE),
    function() seurust::RunUMISampling(umi, 50, FALSE, FALSE)),
  list("RunUMISamplingPerCell", "2,000 × 1,500 counts, 8% nonzero", 5L,
    function() seurat_fn("RunUMISamplingPerCell")(umi, umi_per_cell, FALSE, FALSE),
    function() seurust::RunUMISamplingPerCell(umi, umi_per_cell, FALSE, FALSE)),
  list("ComputeSNN", "1,500 cells, k = 20", 7L,
    function() seurat_fn("ComputeSNN")(nn, 1 / 15),
    function() seurust::ComputeSNN(nn, 1 / 15)),
  list("IntegrateDataC", "800 genes × 200 cells, 60 anchors", 5L,
    function() seurat_fn("IntegrateDataC")(im, weights, expr),
    function() seurust::IntegrateDataC(im, weights, expr)),
  list("FindWeightsC", "400 cells, 600 anchors, k = 10", 5L,
    function() seurat_fn("FindWeightsC")(cells2, distances, anchor_cells2, integration_rownames, cell_index, anchor_score, 0, 1, FALSE),
    function() seurust::FindWeightsC(cells2, distances, anchor_cells2, integration_rownames, cell_index, anchor_score, 0, 1, FALSE)),
  list("ScoreHelper", "300 cells × 16 dimensions", 5L,
    function() seurat_fn("ScoreHelper")(score_snn, query_pca, query_dists, corrected_nns, 8L, FALSE, FALSE),
    function() seurust::ScoreHelper(score_snn, query_pca, query_dists, corrected_nns, 8L, FALSE, FALSE)),
  list("WriteEdgeFile", "400-cell SNN", 3L,
    function() { seurat_fn("WriteEdgeFile")(snn_file, edge_cpp, FALSE); invisible(NULL) },
    function() { seurust::WriteEdgeFile(snn_file, edge_rust, FALSE); invisible(NULL) }),
  list("DirectSNNToFile", "400 cells, k = 15", 3L,
    function() seurat_fn("DirectSNNToFile")(nn_file, 1 / 15, FALSE, direct_cpp),
    function() seurust::DirectSNNToFile(nn_file, 1 / 15, FALSE, direct_rust)),
  list("SNN_SmallestNonzero_Dist", "400 cells × 10 dimensions", 7L,
    function() seurat_fn("SNN_SmallestNonzero_Dist")(snn_file, embed, 5L, nearest_dist),
    function() seurust::SNN_SmallestNonzero_Dist(snn_file, embed, 5L, nearest_dist)),
  list("RunModularityClusteringCpp", "250-cell SNN, 10 iterations", 3L,
    function() seurat_fn("RunModularityClusteringCpp")(mod_snn, 1L, 1, 3L, 1L, 10L, 42L, FALSE, ""),
    function() seurust::RunModularityClusteringCpp(mod_snn, 1L, 1, 3L, 1L, 10L, 42L, FALSE, "")),
  list("fast_dist", "800 × 20 dense, 10 neighbors", 5L,
    function() seurat_fn("fast_dist")(dist_x, dist_y, dist_n),
    function() seurust::fast_dist(dist_x, dist_y, dist_n)),
  list("row_sum_dgcmatrix", "6,000 × 1,500 sparse, 4% nonzero", 11L,
    function() seurat_fn("row_sum_dgcmatrix")(row_x, row_i, row_nr, row_nc),
    function() seurust::row_sum_dgcmatrix(row_x, row_i, row_nr, row_nc)),
  list("row_mean_dgcmatrix", "6,000 × 1,500 sparse, 4% nonzero", 11L,
    function() seurat_fn("row_mean_dgcmatrix")(row_x, row_i, row_nr, row_nc),
    function() seurust::row_mean_dgcmatrix(row_x, row_i, row_nr, row_nc)),
  list("row_var_dgcmatrix", "6,000 × 1,500 sparse, 4% nonzero", 11L,
    function() seurat_fn("row_var_dgcmatrix")(row_x, row_i, row_nr, row_nc),
    function() seurust::row_var_dgcmatrix(row_x, row_i, row_nr, row_nc))
)

rows <- vector("list", length(cases))
cat("| Function | Problem size | Seurat C++ | seurust Rust | Speedup |\n")
cat("| --- | --- | ---: | ---: | ---: |\n")
for (i in seq_along(cases)) {
  case <- cases[[i]]
  name <- case[[1]]
  size <- case[[2]]
  times <- case[[3]]
  message("timing ", name)
  result <- tryCatch(
    time_pair(case[[4]], case[[5]], times = times),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    cat(sprintf("| `%s` | %s | error | error | — |\n", name, size))
    message("  FAILED: ", conditionMessage(result))
    rows[[i]] <- data.frame(
      function_name = name, size = size, cpp_ms = NA, rust_ms = NA,
      speedup = NA, error = conditionMessage(result), stringsAsFactors = FALSE
    )
    next
  }
  speedup_txt <- sprintf("%.2f×", result$speedup)
  cat(sprintf(
    "| `%s` | %s | %s | %s | %s |\n",
    name, size, fmt_ms(result$cpp_us), fmt_ms(result$rust_us), speedup_txt
  ))
  rows[[i]] <- data.frame(
    function_name = name, size = size,
    cpp_ms = result$cpp_us / 1000, rust_ms = result$rust_us / 1000,
    speedup = result$speedup, error = NA_character_, stringsAsFactors = FALSE
  )
  gc(verbose = FALSE)
}

out <- do.call(rbind, rows)
csv_path <- Sys.getenv("SEURUST_BENCH_CSV", unset = file.path(tempdir(), "kernel-bench.csv"))
utils::write.csv(out, csv_path, row.names = FALSE)
cat("\nWrote ", csv_path, "\n", sep = "")
