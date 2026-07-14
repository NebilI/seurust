test_that("seurust loads", {
  expect_true(requireNamespace("seurust", quietly = TRUE))
})

test_that("core kernel helpers are exported", {
  ns <- asNamespace("seurust")
  expect_true(exists("LogNorm", where = ns, mode = "function"))
  expect_true(exists("FastSparseRowScale", where = ns, mode = "function"))
  expect_true(exists("ComputeSNN", where = ns, mode = "function"))
})

test_that("LogNorm returns a dgCMatrix with finite values", {
  skip_if_not_installed("Matrix")
  mat <- Matrix::sparseMatrix(
    i = c(1L, 3L, 2L),
    j = c(1L, 2L, 3L),
    x = c(1, 2, 3),
    dims = c(3L, 3L)
  )
  out <- seurust::LogNorm(mat, scale_factor = 1e4, display_progress = FALSE)
  expect_s4_class(out, "dgCMatrix")
  expect_identical(dim(out), dim(mat))
  expect_true(all(is.finite(out@x)))
  expect_true(all(out@x >= 0))
})

test_that("row_sum_dgcmatrix matches Matrix::rowSums", {
  skip_if_not_installed("Matrix")
  mat <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L, 3L),
    j = c(1L, 2L, 2L, 3L),
    x = c(1, 2, 3, 4),
    dims = c(3L, 3L)
  )
  got <- seurust::row_sum_dgcmatrix(mat@x, mat@i, nrow(mat), ncol(mat))
  expect_equal(as.numeric(got), as.numeric(Matrix::rowSums(mat)), tolerance = 1e-10)
})
