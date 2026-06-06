# Timing comparison: Seurat (C++/Rcpp via Seurat:::) vs SeuratRust.
context("SeuratRust/Seurat timing")

test_that("LogNorm timing", {
  skip_if_no_seuratrust()
  mat <- as(matrix(1:160000, ncol = 400, nrow = 400), "sparseMatrix")
  expect_equal(
    as.matrix(SeuratRust::LogNorm(mat, 1e4, display_progress = FALSE)),
    as.matrix(Seurat:::LogNorm(mat, 1e4, display_progress = FALSE)),
    tolerance = 1e-10
  )
  bench <- benchmark_rust_cpp(
    cpp_fn = function() Seurat:::LogNorm(mat, 1e4, display_progress = FALSE),
    rust_fn = function() SeuratRust::LogNorm(mat, 1e4, display_progress = FALSE),
    n_warmup = 2L,
    n_reps = 10L
  )
  expect_timing_report(bench, "LogNorm")
  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, "LogNorm")
  }
})

test_that("ComputeSNN timing (500 cells)", {
  skip_if_no_seuratrust()
  bench <- benchmark_compute_snn(
    n_cells = 500L,
    n_warmup = 2L,
    n_reps = 10L
  )
  expect_timing_report(bench, attr(bench, "label"))
  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, attr(bench, "label"))
  }
})

test_that("ComputeSNN timing (2000 cells)", {
  skip_if_no_seuratrust()
  bench <- benchmark_compute_snn(
    n_cells = 2000L,
    n_warmup = 1L,
    n_reps = 10L
  )
  expect_timing_report(bench, attr(bench, "label"))
  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, attr(bench, "label"))
  }
})

test_that("row_sum_dgcmatrix timing", {
  skip_if_no_seuratrust()
  mat <- sparseMatrix(
    i = sample.int(3000, 50000, replace = TRUE),
    j = sample.int(800, 50000, replace = TRUE),
    x = runif(50000),
    dims = c(3000L, 800L)
  )
  x <- slot(mat, "x")
  i <- slot(mat, "i")
  nr <- nrow(mat)
  nc <- ncol(mat)
  expect_equal(
    SeuratRust::row_sum_dgcmatrix(x, i, nr, nc),
    Seurat:::row_sum_dgcmatrix(x, i, nr, nc),
    tolerance = 1e-10
  )
  bench <- benchmark_rust_cpp(
    cpp_fn = function() Seurat:::row_sum_dgcmatrix(x, i, nr, nc),
    rust_fn = function() SeuratRust::row_sum_dgcmatrix(x, i, nr, nc),
    n_warmup = 2L,
    n_reps = 10L
  )
  expect_timing_report(bench, "row_sum_dgcmatrix")
  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, "row_sum_dgcmatrix")
  }
})
