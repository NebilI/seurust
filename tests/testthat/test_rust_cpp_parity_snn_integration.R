# Parity tests: SeuratRust vs Seurat (C++/Rcpp) for integration, SNN, and kNN.
# Requires SeuratRust to be installed alongside Seurat.
#
# Run locally:
#   devtools::install("SeuratRust")
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test_rust_cpp_parity_snn_integration.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(testthat)
})

context("SeuratRust/Seurat parity: fast_dist")

test_that("SeuratRust fast_dist matches Seurat fast_dist", {
  skip_if_no_seuratrust()
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
  rust <- SeuratRust::fast_dist(x = x, y = y, n = n)
  expect_equal(cpp, rust, tolerance = 1e-10)
})

context("SeuratRust/Seurat parity: ComputeSNN")

test_that("SeuratRust ComputeSNN matches Seurat ComputeSNN", {
  skip_if_no_seuratrust()
  set.seed(2)
  nn <- matrix(sample(x = 1:6, size = 18, replace = TRUE), nrow = 6, ncol = 3)
  prune <- 0.01
  cpp <- Seurat:::ComputeSNN(nn_ranked = nn, prune = prune)
  rust <- SeuratRust::ComputeSNN(nn_ranked = nn, prune = prune)
  expect_equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10)
})

context("SeuratRust/Seurat parity: IntegrateDataC")

test_that("SeuratRust IntegrateDataC matches Seurat IntegrateDataC", {
  skip_if_no_seuratrust()
  set.seed(3)
  expr <- as(sparseMatrix(
    i = c(0, 1, 2, 0, 1),
    p = c(0, 2, 4, 5),
    x = c(1, 2, 3, 4, 5),
    dims = c(3L, 2L)
  ), "dgCMatrix")
  im <- as(sparseMatrix(
    i = c(0, 1, 0),
    p = c(0, 2, 3),
    x = c(0.5, 0.3, 0.2),
    dims = c(2L, 2L)
  ), "dgCMatrix")
  w <- as(sparseMatrix(
    i = c(0, 1, 0),
    p = c(0, 2, 3),
    x = c(0.4, 0.6, 0.1),
    dims = c(2L, 3L)
  ), "dgCMatrix")
  cpp <- Seurat:::IntegrateDataC(
    integration_matrix = im,
    weights = w,
    expression_cells2 = expr
  )
  rust <- SeuratRust::IntegrateDataC(
    integration_matrix = im,
    weights = w,
    expression_cells2 = expr
  )
  expect_equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10)
})

context("SeuratRust/Seurat parity: FindWeightsC")

test_that("SeuratRust FindWeightsC matches Seurat FindWeightsC (min_dist = 0)", {
  skip_if_no_seuratrust()
  set.seed(4)
  cells2 <- 0:1
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
  rust <- SeuratRust::FindWeightsC(
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

context("SeuratRust/Seurat parity: SNN_SmallestNonzero_Dist")

test_that("SeuratRust SNN_SmallestNonzero_Dist matches Seurat", {
  skip_if_no_seuratrust()
  set.seed(5)
  nn <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  snn <- Seurat:::ComputeSNN(nn_ranked = nn, prune = 0)
  mat <- matrix(rnorm(9), nrow = 3, ncol = 3)
  nearest_dist <- c(0, 0.1, 0)
  cpp <- Seurat:::SNN_SmallestNonzero_Dist(
    snn = snn, mat = mat, n = 2, nearest_dist = nearest_dist
  )
  rust <- SeuratRust::SNN_SmallestNonzero_Dist(
    snn = snn, mat = mat, n = 2, nearest_dist = nearest_dist
  )
  expect_equal(cpp, rust, tolerance = 1e-10)
})

context("SeuratRust/Seurat parity: ScoreHelper")

test_that("SeuratRust ScoreHelper matches Seurat ScoreHelper", {
  skip_if_no_seuratrust()
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
  rust <- SeuratRust::ScoreHelper(
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
