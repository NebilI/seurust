test_that("seurust loads", {
  expect_true(requireNamespace("seurust", quietly = TRUE))
})

test_that("LogNorm is exported", {
  skip_if_not_installed("seurust")
  expect_true(exists("LogNorm", where = asNamespace("seurust"), mode = "function"))
})
