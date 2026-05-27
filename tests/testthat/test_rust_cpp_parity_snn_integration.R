# Parity tests: Rust (_rust) vs C++ (Rcpp) for integration, SNN, and kNN distance code.
# Requires the package to be compiled with Rust/extendr support.
#
# Run locally:
#   devtools::load_all()
#   testthat::test_file("tests/testthat/test_rust_cpp_parity_snn_integration.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(testthat)
})

skip_if_no_rust <- function() {
  if (!exists("compute_snn_rust", mode = "function")) {
    skip("Rust extendr functions not available; compile with Rust toolchain")
  }
}

context("Rust/C++ parity: fast_dist")

test_that("fast_dist_rust matches fast_dist", {
  skip_if_no_rust()
  set.seed(1)
  x <- matrix(rnorm(12), nrow = 4, ncol = 3)
  y <- matrix(rnorm(12), nrow = 4, ncol = 3)
  n <- list(
    c(1, 2, 3),
    c(2, 4, 1),
    c(3, 1, 4),
    c(4, 2, 3)
  )
  cpp <- fast_dist(x = x, y = y, n = n)
  rust <- FastDistRust(x = x, y = y, n = n)
  expect_equal(cpp, rust, tolerance = 1e-10)
})

context("Rust/C++ parity: ComputeSNN")

test_that("ComputeSNNRust matches ComputeSNN", {
  skip_if_no_rust()
  set.seed(2)
  nn <- matrix(sample(x = 1:6, size = 18, replace = TRUE), nrow = 6, ncol = 3)
  prune <- 0.01
  cpp <- ComputeSNN(nn_ranked = nn, prune = prune)
  rust <- ComputeSNNRust(nn_ranked = nn, prune = prune)
  expect_equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10)
})

context("Rust/C++ parity: IntegrateDataC")

test_that("IntegrateDataRust matches IntegrateDataC", {
  skip_if_no_rust()
  set.seed(3)
  # cells x genes layout (matches IntegrateData usage after t(data))
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
  cpp <- IntegrateDataC(
    integration_matrix = im,
    weights = w,
    expression_cells2 = expr
  )
  rust <- IntegrateDataRust(
    integration_matrix = im,
    weights = w,
    expression_cells2 = expr
  )
  expect_equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10)
})

context("Rust/C++ parity: FindWeightsC")

test_that("FindWeightsRust matches FindWeightsC (min_dist = 0)", {
  skip_if_no_rust()
  set.seed(4)
  cells2 <- 0:1
  distances <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, byrow = TRUE)
  anchor_cells2 <- c("a", "b")
  rownames <- c("g1", "g2", "g1")
  cell_index <- matrix(c(1, 2, 2, 1), nrow = 2, byrow = TRUE)
  anchor_score <- c(1, 0.5, 0.8)
  cpp <- FindWeightsC(
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
  rust <- FindWeightsRust(
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

context("Rust/C++ parity: SNN_SmallestNonzero_Dist")

test_that("SNNWidthRust matches SNN_SmallestNonzero_Dist", {
  skip_if_no_rust()
  set.seed(5)
  nn <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  snn <- ComputeSNN(nn_ranked = nn, prune = 0)
  mat <- matrix(rnorm(9), nrow = 3, ncol = 3)
  nearest_dist <- c(0, 0.1, 0)
  cpp <- SNN_SmallestNonzero_Dist(
    snn = snn, mat = mat, n = 2, nearest_dist = nearest_dist
  )
  rust <- SNNWidthRust(
    snn = snn, mat = mat, n = 2, nearest_dist = nearest_dist
  )
  expect_equal(cpp, rust, tolerance = 1e-10)
})

context("Rust/C++ parity: ScoreHelper")

test_that("ScoreHelperRust matches ScoreHelper", {
  skip_if_no_rust()
  set.seed(6)
  nn <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  snn <- ComputeSNN(nn_ranked = nn, prune = 0)
  query_pca <- matrix(rnorm(9), nrow = 3, ncol = 3)
  query_dists <- matrix(abs(rnorm(9)), nrow = 3, ncol = 3)
  corrected_nns <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  cpp <- ScoreHelper(
    snn = snn,
    query_pca = query_pca,
    query_dists = query_dists,
    corrected_nns = corrected_nns,
    k_snn = 2,
    subtract_first_nn = FALSE,
    display_progress = FALSE
  )
  rust <- ScoreHelperRust(
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
