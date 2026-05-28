test_that("SeuratRust loads", {
  expect_true(requireNamespace("SeuratRust", quietly = TRUE))
})

test_that("LogNorm is exported", {
  skip_if_not_installed("SeuratRust")
  expect_true(exists("LogNorm", where = asNamespace("SeuratRust"), mode = "function"))
})
