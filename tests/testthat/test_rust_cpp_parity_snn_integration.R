# Parity tests: RSeurat vs Seurat (C++/Rcpp) for integration, SNN, and kNN.
# Requires RSeurat to be installed alongside Seurat.
#
# Run locally:
#   devtools::install("RSeurat")
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test_rust_cpp_parity_snn_integration.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(testthat)
})

context("RSeurat/Seurat parity: fast_dist")

test_that("RSeurat fast_dist matches Seurat fast_dist", {
  skip_if_no_rseurat()
  set.seed(1)
  x <- matrix(rnorm(12), nrow = 4, ncol = 3)
  y <- matrix(rnorm(12), nrow = 4, ncol = 3)
  n <- list(
    c(1, 2, 3),
    c(2, 4, 1),
    c(3, 1, 4),
    c(4, 2, 3)
  )
  cpp <- Seurat:::fast_dist(x = x, y = y, n = n)
  rust <- RSeurat::fast_dist(x = x, y = y, n = n)
  expect_equal(cpp, rust, tolerance = 1e-10)
})

context("RSeurat/Seurat parity: ComputeSNN")

expect_compute_snn_equal <- function(nn, prune) {
  cpp <- Seurat:::ComputeSNN(nn_ranked = nn, prune = prune)
  rust <- RSeurat::ComputeSNN(nn_ranked = nn, prune = prune)

  expect_s4_class(rust, "dgCMatrix")
  expect_equal(dim(rust), dim(cpp))
  expect_equal(slot(rust, "p"), slot(cpp, "p"))
  expect_equal(slot(rust, "i"), slot(cpp, "i"))
  expect_equal(slot(rust, "x"), slot(cpp, "x"), tolerance = 1e-10)
  expect_equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10)
}

test_that("RSeurat ComputeSNN matches Seurat ComputeSNN", {
  skip_if_no_rseurat()

  no_duplicates <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 1,
      3, 4, 5, 6, 1, 2
    ),
    nrow = 6,
    ncol = 3
  )
  storage.mode(no_duplicates) <- "double"

  duplicate_ranks <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      1, 2, 1, 4, 4, 6,
      2, 3, 3, 5, 5, 1,
      2, 2, 3, 5, 1, 1
    ),
    nrow = 6,
    ncol = 4
  )
  storage.mode(duplicate_ranks) <- "double"

  set.seed(2)
  random_small <- matrix(sample(x = 1:6, size = 18, replace = TRUE), nrow = 6, ncol = 3)
  storage.mode(random_small) <- "double"

  set.seed(22)
  random_larger <- matrix(
    sample.int(200L, 200L * 20L, replace = TRUE),
    nrow = 200L,
    ncol = 20L
  )
  storage.mode(random_larger) <- "double"

  for (nn in list(no_duplicates, duplicate_ranks, random_small, random_larger)) {
    for (prune in c(0, 0.01, 0.5)) {
      expect_compute_snn_equal(nn = nn, prune = prune)
    }
  }
})

context("RSeurat/Seurat parity: IntegrateDataC")

test_that("RSeurat IntegrateDataC matches Seurat IntegrateDataC", {
  skip_if_no_rseurat()
  set.seed(3)
  expr <- as(sparseMatrix(
    i = c(0, 1, 1, 2),
    p = c(0, 2, 4),
    x = c(1, 2, 3, 4),
    dims = c(3L, 2L),
    index1 = FALSE
  ), "dgCMatrix")
  im <- as(sparseMatrix(
    i = c(0, 1, 0),
    p = c(0, 2, 3),
    x = c(0.5, 0.3, 0.2),
    dims = c(2L, 2L),
    index1 = FALSE
  ), "dgCMatrix")
  w <- as(sparseMatrix(
    i = c(0, 1, 0),
    p = c(0, 2, 3),
    x = c(0.4, 0.6, 0.1),
    dims = c(2L, 3L),
    index1 = FALSE
  ), "dgCMatrix")
  cpp <- Seurat:::IntegrateDataC(
    integration_matrix = im,
    weights = w,
    expression_cells2 = expr
  )
  rust <- RSeurat::IntegrateDataC(
    integration_matrix = im,
    weights = w,
    expression_cells2 = expr
  )
  expect_equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10)
})

context("RSeurat/Seurat parity: FindWeightsC")

test_that("RSeurat FindWeightsC matches Seurat FindWeightsC (min_dist = 0)", {
  skip_if_no_rseurat()
  set.seed(4)
  cells2 <- as.numeric(0:1)
  distances <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, byrow = TRUE)
  anchor_cells2 <- c("a", "b")
  rownames <- c("g1", "g2", "g1")
  cell_index <- matrix(c(1, 2, 2, 1), nrow = 2, byrow = TRUE)
  anchor_score <- c(1, 0.5, 0.8)
  cpp <- Seurat:::FindWeightsC(
    cells2 = cells2,
    distances = distances,
    anchor_cells2 = anchor_cells2,
    integration_matrix_rownames = rownames,
    cell_index = cell_index,
    anchor_score = anchor_score,
    min_dist = 0,
    sd = 1,
    display_progress = FALSE
  )
  rust <- RSeurat::FindWeightsC(
    cells2 = cells2,
    distances = distances,
    anchor_cells2 = anchor_cells2,
    integration_matrix_rownames = rownames,
    cell_index = cell_index,
    anchor_score = anchor_score,
    min_dist = 0,
    sd = 1,
    display_progress = FALSE
  )
  expect_equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10)
})

context("RSeurat/Seurat parity: SNN_SmallestNonzero_Dist")

test_that("RSeurat SNN_SmallestNonzero_Dist matches Seurat", {
  skip_if_no_rseurat()
  set.seed(5)
  nn <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  snn <- Seurat:::ComputeSNN(nn_ranked = nn, prune = 0)
  mat <- matrix(rnorm(9), nrow = 3, ncol = 3)
  nearest_dist <- c(0, 0.1, 0)
  cpp <- Seurat:::SNN_SmallestNonzero_Dist(
    snn = snn, mat = mat, n = 2, nearest_dist = nearest_dist
  )
  rust <- RSeurat::SNN_SmallestNonzero_Dist(
    snn = snn, mat = mat, n = 2, nearest_dist = nearest_dist
  )
  expect_equal(cpp, rust, tolerance = 1e-10)
})

context("RSeurat/Seurat parity: ScoreHelper")

test_that("RSeurat ScoreHelper matches Seurat ScoreHelper", {
  skip_if_no_rseurat()
  set.seed(6)
  nn <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  snn <- Seurat:::ComputeSNN(nn_ranked = nn, prune = 0)
  query_pca <- matrix(rnorm(9), nrow = 3, ncol = 3)
  query_dists <- matrix(abs(rnorm(9)), nrow = 3, ncol = 3)
  corrected_nns <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  cpp <- Seurat:::ScoreHelper(
    snn = snn,
    query_pca = query_pca,
    query_dists = query_dists,
    corrected_nns = corrected_nns,
    k_snn = 2,
    subtract_first_nn = FALSE,
    display_progress = FALSE
  )
  rust <- RSeurat::ScoreHelper(
    snn = snn,
    query_pca = query_pca,
    query_dists = query_dists,
    corrected_nns = corrected_nns,
    k_snn = 2,
    subtract_first_nn = FALSE,
    display_progress = FALSE
  )
  expect_equal(cpp, rust, tolerance = 1e-10)
})
