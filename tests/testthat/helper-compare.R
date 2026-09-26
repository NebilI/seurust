#' Shared fixtures and helpers for seurust vs Seurat (C++) parity and comparison.

SEURAT_NATIVE_FUNCTIONS <- c(
  "ComputeSNN",
  "DirectSNNToFile",
  "FastCov",
  "FastCovMats",
  "FastExpMean",
  "FastLogVMR",
  "FastRBind",
  "FastSparseRowScale",
  "FastSparseRowScaleWithKnownStats",
  "FindWeightsC",
  "GraphToNeighborHelper",
  "IntegrateDataC",
  "LogNorm",
  "ReplaceColsC",
  "RowMergeMatrices",
  "RowVar",
  "RunModularityClusteringCpp",
  "RunUMISampling",
  "RunUMISamplingPerCell",
  "SNN_SmallestNonzero_Dist",
  "ScoreHelper",
  "SparseRowVar",
  "SparseRowVar2",
  "SparseRowVarStd",
  "Standardize",
  "WriteEdgeFile",
  "fast_dist",
  "row_mean_dgcmatrix",
  "row_sum_dgcmatrix",
  "row_var_dgcmatrix"
)

compare_n_reps <- function(default = 10L) {
  env <- Sys.getenv("SEURUST_BENCH_REPS", unset = "")
  if (nzchar(env)) as.integer(env) else as.integer(default)
}

compare_n_mem <- function(default = 4L) {
  env <- Sys.getenv("SEURUST_BENCH_MEM_REPS", unset = "")
  if (nzchar(env)) as.integer(env) else as.integer(default)
}

as_dgc <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  if (is.list(x) && !is.null(x$x) && !is.null(x$i) && !is.null(x$p) && !is.null(x$Dim)) {
    return(methods::new(
      Class = "dgCMatrix",
      x = as.numeric(x$x),
      i = as.integer(x$i),
      p = as.integer(x$p),
      Dim = as.integer(x$Dim)
    ))
  }
  as(x, "dgCMatrix")
}

expect_numeric_equal <- function(a, b, tolerance = 1e-10) {
  testthat::expect_equal(as.numeric(a), as.numeric(b), tolerance = tolerance)
}

expect_matrix_equal <- function(a, b, tolerance = 1e-10) {
  testthat::expect_equal(as.matrix(a), as.matrix(b), tolerance = tolerance, check.attributes = FALSE)
}

expect_dgc_equal <- function(a, b, tolerance = 1e-10) {
  expect_matrix_equal(as_dgc(a), as_dgc(b), tolerance = tolerance)
}

expect_list_matrices_equal <- function(a, b, tolerance = 1e-10) {
  testthat::expect_equal(length(a), length(b))
  for (i in seq_along(a)) {
    expect_matrix_equal(a[[i]], b[[i]], tolerance = tolerance)
  }
}

read_edge_table <- function(path) {
  if (!file.exists(path) || is.na(file.size(path)) || file.size(path) == 0) {
    return(data.frame(col = integer(), row = integer(), value = numeric()))
  }
  tab <- utils::read.table(
    path,
    header = FALSE,
    col.names = c("col", "row", "value"),
    stringsAsFactors = FALSE
  )
  tab[order(tab$col, tab$row), , drop = FALSE]
}

expect_edge_files_equal <- function(path_a, path_b, tolerance = 1e-10) {
  a <- read_edge_table(path_a)
  b <- read_edge_table(path_b)
  testthat::expect_equal(nrow(a), nrow(b))
  testthat::expect_equal(a$col, b$col)
  testthat::expect_equal(a$row, b$row)
  testthat::expect_equal(a$value, b$value, tolerance = tolerance)
}

stop_if_not_equal <- function(a, b, tolerance = 1e-10, label = "values") {
  ok <- isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
  if (!ok) {
    stop(label, " do not match between Seurat and seurust", call. = FALSE)
  }
  invisible(TRUE)
}

make_sparse <- function(nrow, ncol, density = 0.12, seed = 1L) {
  set.seed(seed)
  as(Matrix::rsparsematrix(nrow, ncol, density = density, rand.x = stats::runif), "dgCMatrix")
}

make_dense <- function(nrow, ncol, seed = 1L) {
  set.seed(seed)
  matrix(stats::rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
}

make_count_sparse <- function(nrow, ncol, density = 0.2, seed = 1L) {
  set.seed(seed)
  as(
    Matrix::rsparsematrix(
      nrow,
      ncol,
      density = density,
      rand.x = function(n) stats::rpois(n, lambda = 4) + 1
    ),
    "dgCMatrix"
  )
}

# Symmetric regular ring graph (same nnz per row/column). Seurat's
# GraphToNeighborHelper is used on symmetric graphs such as SNN.
make_regular_knn_graph <- function(n = 8L, k = 3L, seed = 1L) {
  set.seed(seed)
  stopifnot(k >= 1L, n > 2L * k)
  i <- integer(0)
  j <- integer(0)
  x <- numeric(0)
  for (row in seq_len(n)) {
    offs <- c(-seq_len(k), seq_len(k))
    nbrs <- ((row - 1L + offs) %% n) + 1L
    i <- c(i, rep(row, length(nbrs)))
    j <- c(j, nbrs)
    # Symmetric weights so (i, j) == (j, i); unique per pair for stable sorts.
    pair <- pmin(row, nbrs) * (n + 1L) + pmax(row, nbrs)
    x <- c(x, abs(as.numeric(offs)) + pair / 10000)
  }
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
}

modularity_toy_snn <- function() {
  node1 <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
    1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 8, 8, 8, 9, 13,
    14, 14, 15, 15, 18, 18, 19, 20, 20, 22, 22, 23, 23, 23, 23, 23, 24,
    24, 24, 25, 26, 26, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32
  )
  node2 <- c(
    1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 17, 19, 21, 31, 2, 3, 7, 13,
    17, 19, 21, 30, 3, 7, 8, 9, 13, 27, 28, 32, 7, 12, 13, 6, 10, 6, 10,
    16, 16, 30, 32, 33, 33, 33, 32, 33, 32, 33, 32, 33, 33, 32, 33, 32,
    33, 25, 27, 29, 32, 33, 25, 27, 31, 31, 29, 33, 33, 31, 33, 32, 33,
    32, 33, 32, 33, 33
  )
  Matrix::sparseMatrix(i = node2 + 1, j = node1 + 1, x = 1.0)
}

kernel_case <- function(name, size, cpp_fn, rust_fn, check_parity) {
  list(
    name = name,
    size = size,
    cpp_fn = cpp_fn,
    rust_fn = rust_fn,
    check_parity = check_parity
  )
}

#' Catalog of every ported native kernel, with Seurat and seurust callables.
#' @keywords internal
seurust_kernel_cases <- function() {
  log_mat <- make_sparse(200L, 200L, density = 0.15, seed = 11L)
  dense_mat <- make_dense(300L, 60L, seed = 12L)
  cov_mat <- dense_mat[, seq_len(40L), drop = FALSE]
  cov_mat2 <- make_dense(300L, 30L, seed = 13L)
  rbind_a <- dense_mat[seq_len(80L), seq_len(20L), drop = FALSE]
  rbind_b <- dense_mat[81:140, seq_len(20L), drop = FALSE]
  sparse_stats <- make_sparse(250L, 180L, density = 0.14, seed = 14L)
  mu <- Matrix::rowMeans(sparse_stats)
  sigma <- sqrt(pmax(apply(as.matrix(sparse_stats), 1L, stats::var), .Machine$double.eps))
  scale_mat <- make_sparse(180L, 220L, density = 0.12, seed = 15L)
  scale_mu <- Matrix::rowMeans(scale_mat)
  scale_sd <- sqrt(pmax(apply(as.matrix(scale_mat), 1L, stats::var), .Machine$double.eps))
  merge_a <- as(make_sparse(40L, 30L, density = 0.2, seed = 16L), "RsparseMatrix")
  merge_b <- as(make_sparse(40L, 30L, density = 0.2, seed = 17L), "RsparseMatrix")
  merge_n1 <- paste0("g", sample.int(60L, 40L))
  merge_n2 <- paste0("g", sample.int(60L, 40L))
  merge_all <- union(merge_n1, merge_n2)
  replace_mat <- make_sparse(30L, 80L, density = 0.15, seed = 18L)
  replace_src <- make_sparse(30L, 8L, density = 0.2, seed = 19L)
  replace_idx <- as.numeric(c(0L, 3L, 7L, 10L, 22L, 40L, 55L, 70L))
  knn_graph <- make_regular_knn_graph(n = 24L, k = 2L, seed = 20L)
  umi_mat <- make_count_sparse(80L, 60L, density = 0.25, seed = 21L)
  umi_per_cell <- as.numeric(rep(20, ncol(umi_mat)))
  nn_small <- make_compute_snn_nn(n_cells = 120L, k = 12L, seed = 22L)
  snn_small <- Seurat:::ComputeSNN(nn_small, prune = 0.01)
  # IntegrateDataC computes expression - t(weights) %*% integration_matrix,
  # so weights is K x G, integration_matrix is K x C2, expression is G x C2.
  expr <- make_sparse(40L, 25L, density = 0.2, seed = 23L)
  im <- make_sparse(30L, 25L, density = 0.15, seed = 24L)
  weights <- make_sparse(30L, 40L, density = 0.1, seed = 25L)
  n_cells_w <- 80L
  n_anchors <- 120L
  k_weight <- 8L
  cells2 <- as.numeric(seq_len(n_cells_w) - 1L)
  distances <- matrix(stats::runif(n_cells_w * k_weight), nrow = n_cells_w)
  cell_index <- matrix(
    sample.int(n_anchors, n_cells_w * k_weight, replace = TRUE),
    nrow = n_cells_w
  )
  storage.mode(cell_index) <- "double"
  anchor_cells2 <- paste0("cell", sample.int(40L, n_anchors, replace = TRUE))
  integration_rownames <- paste0("cell", sample.int(40L, n_anchors, replace = TRUE))
  anchor_score <- stats::runif(n_anchors, min = 0.1, max = 1)
  score_nn <- make_compute_snn_nn(n_cells = 60L, k = 8L, seed = 26L)
  score_snn <- Seurat:::ComputeSNN(score_nn, prune = 0)
  query_pca <- make_dense(12L, 60L, seed = 27L)
  query_dists <- abs(make_dense(60L, 8L, seed = 28L))
  corrected_nns <- make_compute_snn_nn(n_cells = 60L, k = 8L, seed = 29L)
  dist_x <- make_dense(80L, 12L, seed = 30L)
  dist_y <- make_dense(80L, 12L, seed = 31L)
  dist_n <- replicate(80L, as.numeric(sample.int(80L, 6L)), simplify = FALSE)
  row_mat <- make_sparse(400L, 120L, density = 0.08, seed = 32L)
  row_x <- methods::slot(row_mat, "x")
  row_i <- methods::slot(row_mat, "i")
  row_nr <- nrow(row_mat)
  row_nc <- ncol(row_mat)
  snn_mod <- modularity_toy_snn()
  modularity_args <- list(
    modularityFunction = 1L,
    resolution = 1.0,
    algorithm = 3L,
    nRandomStarts = 2L,
    nIterations = 10L,
    randomSeed = 42L,
    printOutput = FALSE,
    edgefilename = ""
  )
  nearest_dist <- abs(stats::rnorm(nrow(snn_small)))
  nearest_dist[c(1L, 5L)] <- 0
  embed <- make_dense(nrow(snn_small), 8L, seed = 33L)

  list(
    kernel_case(
      "LogNorm",
      "200x200 sparse",
      function() Seurat:::LogNorm(log_mat, 1e4, display_progress = FALSE),
      function() seurust::LogNorm(log_mat, 1e4, display_progress = FALSE),
      function() expect_dgc_equal(
        Seurat:::LogNorm(log_mat, 1e4, display_progress = FALSE),
        seurust::LogNorm(log_mat, 1e4, display_progress = FALSE)
      )
    ),
    kernel_case(
      "Standardize",
      "300x60 dense",
      function() Seurat:::Standardize(dense_mat, display_progress = FALSE),
      function() seurust::Standardize(dense_mat, display_progress = FALSE),
      function() expect_matrix_equal(
        Seurat:::Standardize(dense_mat, display_progress = FALSE),
        seurust::Standardize(dense_mat, display_progress = FALSE)
      )
    ),
    kernel_case(
      "FastCov",
      "300x40 dense",
      function() Seurat:::FastCov(cov_mat, center = TRUE),
      function() seurust::FastCov(cov_mat, center = TRUE),
      function() expect_matrix_equal(
        Seurat:::FastCov(cov_mat, center = TRUE),
        seurust::FastCov(cov_mat, center = TRUE)
      )
    ),
    kernel_case(
      "FastCovMats",
      "300x40 vs 300x30",
      function() Seurat:::FastCovMats(cov_mat, cov_mat2, center = TRUE),
      function() seurust::FastCovMats(cov_mat, cov_mat2, center = TRUE),
      function() expect_matrix_equal(
        Seurat:::FastCovMats(cov_mat, cov_mat2, center = TRUE),
        seurust::FastCovMats(cov_mat, cov_mat2, center = TRUE)
      )
    ),
    kernel_case(
      "FastRBind",
      "80x20 + 60x20",
      function() Seurat:::FastRBind(rbind_a, rbind_b),
      function() seurust::FastRBind(rbind_a, rbind_b),
      function() expect_matrix_equal(
        Seurat:::FastRBind(rbind_a, rbind_b),
        seurust::FastRBind(rbind_a, rbind_b)
      )
    ),
    kernel_case(
      "RowVar",
      "300x60 dense",
      function() Seurat:::RowVar(dense_mat),
      function() seurust::RowVar(dense_mat),
      function() expect_numeric_equal(
        Seurat:::RowVar(dense_mat),
        seurust::RowVar(dense_mat)
      )
    ),
    kernel_case(
      "FastExpMean",
      "250x180 sparse",
      function() Seurat:::FastExpMean(sparse_stats, display_progress = FALSE),
      function() seurust::FastExpMean(sparse_stats, display_progress = FALSE),
      function() expect_numeric_equal(
        Seurat:::FastExpMean(sparse_stats, display_progress = FALSE),
        seurust::FastExpMean(sparse_stats, display_progress = FALSE)
      )
    ),
    kernel_case(
      "SparseRowVar",
      "250x180 sparse",
      function() Seurat:::SparseRowVar(sparse_stats, display_progress = FALSE),
      function() seurust::SparseRowVar(sparse_stats, display_progress = FALSE),
      function() expect_numeric_equal(
        Seurat:::SparseRowVar(sparse_stats, display_progress = FALSE),
        seurust::SparseRowVar(sparse_stats, display_progress = FALSE)
      )
    ),
    kernel_case(
      "SparseRowVar2",
      "250x180 sparse",
      function() Seurat:::SparseRowVar2(sparse_stats, mu = mu, display_progress = FALSE),
      function() seurust::SparseRowVar2(sparse_stats, mu = mu, display_progress = FALSE),
      function() expect_numeric_equal(
        Seurat:::SparseRowVar2(sparse_stats, mu = mu, display_progress = FALSE),
        seurust::SparseRowVar2(sparse_stats, mu = mu, display_progress = FALSE)
      )
    ),
    kernel_case(
      "SparseRowVarStd",
      "250x180 sparse",
      function() {
        Seurat:::SparseRowVarStd(
          sparse_stats, mu = mu, sd = sigma, vmax = 10, display_progress = FALSE
        )
      },
      function() {
        seurust::SparseRowVarStd(
          sparse_stats, mu = mu, sd = sigma, vmax = 10, display_progress = FALSE
        )
      },
      function() expect_numeric_equal(
        Seurat:::SparseRowVarStd(
          sparse_stats, mu = mu, sd = sigma, vmax = 10, display_progress = FALSE
        ),
        seurust::SparseRowVarStd(
          sparse_stats, mu = mu, sd = sigma, vmax = 10, display_progress = FALSE
        )
      )
    ),
    kernel_case(
      "FastLogVMR",
      "250x180 sparse",
      function() Seurat:::FastLogVMR(sparse_stats, display_progress = FALSE),
      function() seurust::FastLogVMR(sparse_stats, display_progress = FALSE),
      function() expect_numeric_equal(
        Seurat:::FastLogVMR(sparse_stats, display_progress = FALSE),
        seurust::FastLogVMR(sparse_stats, display_progress = FALSE)
      )
    ),
    kernel_case(
      "FastSparseRowScale",
      "180x220 sparse",
      function() {
        Seurat:::FastSparseRowScale(
          scale_mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        )
      },
      function() {
        seurust::FastSparseRowScale(
          scale_mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        )
      },
      function() expect_matrix_equal(
        Seurat:::FastSparseRowScale(
          scale_mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        ),
        seurust::FastSparseRowScale(
          scale_mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        )
      )
    ),
    kernel_case(
      "FastSparseRowScaleWithKnownStats",
      "180x220 sparse",
      function() {
        Seurat:::FastSparseRowScaleWithKnownStats(
          scale_mat, mu = scale_mu, sigma = scale_sd,
          scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        )
      },
      function() {
        seurust::FastSparseRowScaleWithKnownStats(
          scale_mat, mu = scale_mu, sigma = scale_sd,
          scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        )
      },
      function() expect_matrix_equal(
        Seurat:::FastSparseRowScaleWithKnownStats(
          scale_mat, mu = scale_mu, sigma = scale_sd,
          scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        ),
        seurust::FastSparseRowScaleWithKnownStats(
          scale_mat, mu = scale_mu, sigma = scale_sd,
          scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE
        )
      )
    ),
    kernel_case(
      "RowMergeMatrices",
      "40x30 + 40x30 CSR",
      function() {
        Seurat:::RowMergeMatrices(
          mat1 = merge_a, mat2 = merge_b,
          mat1_rownames = merge_n1, mat2_rownames = merge_n2,
          all_rownames = merge_all
        )
      },
      function() {
        seurust::RowMergeMatrices(
          mat1 = merge_a, mat2 = merge_b,
          mat1_rownames = merge_n1, mat2_rownames = merge_n2,
          all_rownames = merge_all
        )
      },
      function() expect_dgc_equal(
        Seurat:::RowMergeMatrices(
          mat1 = merge_a, mat2 = merge_b,
          mat1_rownames = merge_n1, mat2_rownames = merge_n2,
          all_rownames = merge_all
        ),
        seurust::RowMergeMatrices(
          mat1 = merge_a, mat2 = merge_b,
          mat1_rownames = merge_n1, mat2_rownames = merge_n2,
          all_rownames = merge_all
        )
      )
    ),
    kernel_case(
      "ReplaceColsC",
      "30x80 replace 8 cols",
      function() Seurat:::ReplaceColsC(replace_mat, replace_idx, replace_src),
      function() seurust::ReplaceColsC(replace_mat, replace_idx, replace_src),
      function() expect_dgc_equal(
        Seurat:::ReplaceColsC(replace_mat, replace_idx, replace_src),
        seurust::ReplaceColsC(replace_mat, replace_idx, replace_src)
      )
    ),
    kernel_case(
      "GraphToNeighborHelper",
      "24 cells, k=2 (symmetric)",
      function() Seurat:::GraphToNeighborHelper(knn_graph),
      function() seurust::GraphToNeighborHelper(knn_graph),
      function() expect_list_matrices_equal(
        Seurat:::GraphToNeighborHelper(knn_graph),
        seurust::GraphToNeighborHelper(knn_graph)
      )
    ),
    kernel_case(
      "RunUMISampling",
      "80x60 counts",
      function() {
        set.seed(101)
        Seurat:::RunUMISampling(umi_mat, sample_val = 25L, upsample = FALSE, display_progress = FALSE)
      },
      function() {
        set.seed(101)
        seurust::RunUMISampling(umi_mat, sample_val = 25L, upsample = FALSE, display_progress = FALSE)
      },
      function() {
        set.seed(101)
        cpp <- Seurat:::RunUMISampling(umi_mat, sample_val = 25L, upsample = FALSE, display_progress = FALSE)
        set.seed(101)
        rust <- seurust::RunUMISampling(umi_mat, sample_val = 25L, upsample = FALSE, display_progress = FALSE)
        expect_dgc_equal(cpp, rust)
      }
    ),
    kernel_case(
      "RunUMISamplingPerCell",
      "80x60 counts",
      function() {
        set.seed(102)
        Seurat:::RunUMISamplingPerCell(
          umi_mat, sample_val = umi_per_cell, upsample = FALSE, display_progress = FALSE
        )
      },
      function() {
        set.seed(102)
        seurust::RunUMISamplingPerCell(
          umi_mat, sample_val = umi_per_cell, upsample = FALSE, display_progress = FALSE
        )
      },
      function() {
        set.seed(102)
        cpp <- Seurat:::RunUMISamplingPerCell(
          umi_mat, sample_val = umi_per_cell, upsample = FALSE, display_progress = FALSE
        )
        set.seed(102)
        rust <- seurust::RunUMISamplingPerCell(
          umi_mat, sample_val = umi_per_cell, upsample = FALSE, display_progress = FALSE
        )
        expect_dgc_equal(cpp, rust)
      }
    ),
    kernel_case(
      "ComputeSNN",
      "120 cells, k=12",
      function() Seurat:::ComputeSNN(nn_small, prune = 0.01),
      function() seurust::ComputeSNN(nn_small, prune = 0.01),
      function() expect_dgc_equal(
        Seurat:::ComputeSNN(nn_small, prune = 0.01),
        seurust::ComputeSNN(nn_small, prune = 0.01)
      )
    ),
    kernel_case(
      "IntegrateDataC",
      "40 genes, 25 cells, 30 anchors",
      function() Seurat:::IntegrateDataC(im, weights, expr),
      function() seurust::IntegrateDataC(im, weights, expr),
      function() expect_dgc_equal(
        Seurat:::IntegrateDataC(im, weights, expr),
        seurust::IntegrateDataC(im, weights, expr)
      )
    ),
    kernel_case(
      "FindWeightsC",
      "80 cells, 120 anchors",
      function() {
        Seurat:::FindWeightsC(
          cells2 = cells2, distances = distances, anchor_cells2 = anchor_cells2,
          integration_matrix_rownames = integration_rownames, cell_index = cell_index,
          anchor_score = anchor_score, min_dist = 0, sd = 1, display_progress = FALSE
        )
      },
      function() {
        seurust::FindWeightsC(
          cells2 = cells2, distances = distances, anchor_cells2 = anchor_cells2,
          integration_matrix_rownames = integration_rownames, cell_index = cell_index,
          anchor_score = anchor_score, min_dist = 0, sd = 1, display_progress = FALSE
        )
      },
      function() expect_dgc_equal(
        Seurat:::FindWeightsC(
          cells2 = cells2, distances = distances, anchor_cells2 = anchor_cells2,
          integration_matrix_rownames = integration_rownames, cell_index = cell_index,
          anchor_score = anchor_score, min_dist = 0, sd = 1, display_progress = FALSE
        ),
        seurust::FindWeightsC(
          cells2 = cells2, distances = distances, anchor_cells2 = anchor_cells2,
          integration_matrix_rownames = integration_rownames, cell_index = cell_index,
          anchor_score = anchor_score, min_dist = 0, sd = 1, display_progress = FALSE
        )
      )
    ),
    kernel_case(
      "ScoreHelper",
      "60 cells x 12 dims",
      function() {
        Seurat:::ScoreHelper(
          snn = score_snn, query_pca = query_pca, query_dists = query_dists,
          corrected_nns = corrected_nns, k_snn = 8L, subtract_first_nn = FALSE,
          display_progress = FALSE
        )
      },
      function() {
        seurust::ScoreHelper(
          snn = score_snn, query_pca = query_pca, query_dists = query_dists,
          corrected_nns = corrected_nns, k_snn = 8L, subtract_first_nn = FALSE,
          display_progress = FALSE
        )
      },
      function() expect_numeric_equal(
        Seurat:::ScoreHelper(
          snn = score_snn, query_pca = query_pca, query_dists = query_dists,
          corrected_nns = corrected_nns, k_snn = 8L, subtract_first_nn = FALSE,
          display_progress = FALSE
        ),
        seurust::ScoreHelper(
          snn = score_snn, query_pca = query_pca, query_dists = query_dists,
          corrected_nns = corrected_nns, k_snn = 8L, subtract_first_nn = FALSE,
          display_progress = FALSE
        )
      )
    ),
    kernel_case(
      "WriteEdgeFile",
      "120-cell SNN",
      function() {
        path <- tempfile(fileext = ".edgelist")
        on.exit(unlink(path), add = TRUE)
        Seurat:::WriteEdgeFile(snn_small, path, display_progress = FALSE)
        file.exists(path)
      },
      function() {
        path <- tempfile(fileext = ".edgelist")
        on.exit(unlink(path), add = TRUE)
        seurust::WriteEdgeFile(snn_small, path, display_progress = FALSE)
        file.exists(path)
      },
      function() {
        cpp_path <- tempfile(fileext = ".edgelist")
        rust_path <- tempfile(fileext = ".edgelist")
        on.exit(unlink(c(cpp_path, rust_path)), add = TRUE)
        Seurat:::WriteEdgeFile(snn_small, cpp_path, display_progress = FALSE)
        seurust::WriteEdgeFile(snn_small, rust_path, display_progress = FALSE)
        expect_edge_files_equal(cpp_path, rust_path)
      }
    ),
    kernel_case(
      "DirectSNNToFile",
      "120 cells, k=12",
      function() {
        path <- tempfile(fileext = ".edgelist")
        on.exit(unlink(path), add = TRUE)
        Seurat:::DirectSNNToFile(nn_small, 0.01, FALSE, path)
      },
      function() {
        path <- tempfile(fileext = ".edgelist")
        on.exit(unlink(path), add = TRUE)
        seurust::DirectSNNToFile(nn_small, 0.01, FALSE, path)
      },
      function() {
        cpp_path <- tempfile(fileext = ".edgelist")
        rust_path <- tempfile(fileext = ".edgelist")
        on.exit(unlink(c(cpp_path, rust_path)), add = TRUE)
        cpp <- Seurat:::DirectSNNToFile(nn_small, 0.01, FALSE, cpp_path)
        rust <- seurust::DirectSNNToFile(nn_small, 0.01, FALSE, rust_path)
        expect_dgc_equal(cpp, rust)
        expect_edge_files_equal(cpp_path, rust_path)
      }
    ),
    kernel_case(
      "SNN_SmallestNonzero_Dist",
      "120 cells x 8 dims",
      function() {
        Seurat:::SNN_SmallestNonzero_Dist(
          snn = snn_small, mat = embed, n = 5L, nearest_dist = nearest_dist
        )
      },
      function() {
        seurust::SNN_SmallestNonzero_Dist(
          snn = snn_small, mat = embed, n = 5L, nearest_dist = nearest_dist
        )
      },
      function() expect_numeric_equal(
        Seurat:::SNN_SmallestNonzero_Dist(
          snn = snn_small, mat = embed, n = 5L, nearest_dist = nearest_dist
        ),
        seurust::SNN_SmallestNonzero_Dist(
          snn = snn_small, mat = embed, n = 5L, nearest_dist = nearest_dist
        )
      )
    ),
    kernel_case(
      "RunModularityClusteringCpp",
      "34-node SNN",
      function() {
        do.call(Seurat:::RunModularityClusteringCpp, c(list(SNN = snn_mod), modularity_args))
      },
      function() {
        do.call(seurust::RunModularityClusteringCpp, c(list(SNN = snn_mod), modularity_args))
      },
      function() {
        cpp <- do.call(Seurat:::RunModularityClusteringCpp, c(list(SNN = snn_mod), modularity_args))
        rust <- do.call(seurust::RunModularityClusteringCpp, c(list(SNN = snn_mod), modularity_args))
        testthat::expect_equal(rust, cpp)
      }
    ),
    kernel_case(
      "fast_dist",
      "80 cells x 12 dims, k=6",
      function() Seurat:::fast_dist(x = dist_x, y = dist_y, n = dist_n),
      function() seurust::fast_dist(x = dist_x, y = dist_y, n = dist_n),
      function() {
        testthat::expect_equal(
          Seurat:::fast_dist(x = dist_x, y = dist_y, n = dist_n),
          seurust::fast_dist(x = dist_x, y = dist_y, n = dist_n),
          tolerance = 1e-10
        )
      }
    ),
    kernel_case(
      "row_sum_dgcmatrix",
      "400x120 sparse",
      function() Seurat:::row_sum_dgcmatrix(row_x, row_i, row_nr, row_nc),
      function() seurust::row_sum_dgcmatrix(row_x, row_i, row_nr, row_nc),
      function() expect_numeric_equal(
        Seurat:::row_sum_dgcmatrix(row_x, row_i, row_nr, row_nc),
        seurust::row_sum_dgcmatrix(row_x, row_i, row_nr, row_nc)
      )
    ),
    kernel_case(
      "row_mean_dgcmatrix",
      "400x120 sparse",
      function() Seurat:::row_mean_dgcmatrix(row_x, row_i, row_nr, row_nc),
      function() seurust::row_mean_dgcmatrix(row_x, row_i, row_nr, row_nc),
      function() expect_numeric_equal(
        Seurat:::row_mean_dgcmatrix(row_x, row_i, row_nr, row_nc),
        seurust::row_mean_dgcmatrix(row_x, row_i, row_nr, row_nc)
      )
    ),
    kernel_case(
      "row_var_dgcmatrix",
      "400x120 sparse",
      function() Seurat:::row_var_dgcmatrix(row_x, row_i, row_nr, row_nc),
      function() seurust::row_var_dgcmatrix(row_x, row_i, row_nr, row_nc),
      function() expect_numeric_equal(
        Seurat:::row_var_dgcmatrix(row_x, row_i, row_nr, row_nc),
        seurust::row_var_dgcmatrix(row_x, row_i, row_nr, row_nc)
      )
    )
  )
}

run_kernel_parity <- function(case) {
  case$check_parity()
  invisible(TRUE)
}

run_kernel_compare <- function(case, n_warmup = 1L, n_reps = compare_n_reps(), n_mem = compare_n_mem()) {
  parity_ok <- TRUE
  parity_error <- NULL
  tryCatch(
    case$check_parity(),
    error = function(e) {
      parity_ok <<- FALSE
      parity_error <<- conditionMessage(e)
    }
  )
  cmp <- compare_rust_cpp(
    cpp_fn = case$cpp_fn,
    rust_fn = case$rust_fn,
    n_warmup = n_warmup,
    n_reps = n_reps,
    n_mem = n_mem
  )
  register_compare_row(
    name = case$name,
    size = case$size,
    parity = parity_ok,
    cmp = cmp
  )
  list(parity_ok = parity_ok, parity_error = parity_error, cmp = cmp)
}
