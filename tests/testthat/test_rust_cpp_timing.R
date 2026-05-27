# Timing comparison for other C++ vs Rust ports (pure Rust implementations).
context("Rust/C++ timing")

test_that("LogNorm timing", {
  skip_on_cran()
  mat <- as(matrix(1:160000, ncol = 400, nrow = 400), "sparseMatrix")
  cpp_fn <- function() LogNorm(mat, 1e4, display_progress = FALSE)
  rust_fn <- function() LogNormRust(mat, 1e4, display_progress = FALSE)
  stopifnot(all.equal(as.matrix(cpp_fn()), as.matrix(rust_fn()), tolerance = 1e-10))

  bench <- benchmark_rust_cpp(cpp_fn, rust_fn, n_warmup = 2L, n_reps = 10L)
  message(format_benchmark(bench, "LogNorm (400x400 sparse)"))

  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, "LogNorm")
  }
})

test_that("ComputeSNN timing", {
  skip_on_cran()
  set.seed(1)
  nn <- matrix(
    sample.int(500, 500 * 20, replace = TRUE),
    nrow = 500,
    ncol = 20
  )
  storage.mode(nn) <- "double"
  cpp_fn <- function() ComputeSNN(nn, 0.01)
  rust_fn <- function() Seurat:::ComputeSNNRust(nn, 0.01)
  stopifnot(all.equal(as.matrix(cpp_fn()), as.matrix(rust_fn()), tolerance = 1e-10))

  bench <- benchmark_rust_cpp(cpp_fn, rust_fn, n_warmup = 2L, n_reps = 10L)
  message(format_benchmark(bench, "ComputeSNN (500 cells, k=20)"))

  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, "ComputeSNN")
  }
})

test_that("Sparse row sum timing", {
  skip_on_cran()
  mat <- Matrix::sparseMatrix(
    i = sample.int(3000, 50000, replace = TRUE),
    j = sample.int(800, 50000, replace = TRUE),
    x = runif(50000),
    dims = c(3000L, 800L)
  )
  x <- slot(mat, "x")
  i <- slot(mat, "i")
  cpp_fn <- function() Seurat:::row_sum_dgcmatrix(x, i, nrow(mat), ncol(mat))
  rust_fn <- function() Seurat:::row_sum_dgcmatrix_rust(x, i, nrow(mat), ncol(mat))
  stopifnot(all.equal(cpp_fn(), rust_fn()))

  bench <- benchmark_rust_cpp(cpp_fn, rust_fn, n_warmup = 2L, n_reps = 10L)
  message(format_benchmark(bench, "row_sum_dgcmatrix (3000x800 sparse)"))

  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, "row_sum_dgcmatrix")
  }
})
