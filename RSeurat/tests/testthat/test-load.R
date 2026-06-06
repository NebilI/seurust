test_that("RSeurat loads", {
  expect_true(requireNamespace("RSeurat", quietly = TRUE))
})

test_that("LogNorm is exported", {
  skip_if_not_installed("RSeurat")
  expect_true(exists("LogNorm", where = asNamespace("RSeurat"), mode = "function"))
})
