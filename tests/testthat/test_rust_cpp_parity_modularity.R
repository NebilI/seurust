# C++ (Seurat) vs Rust (seurust) parity for modularity clustering.
context("ModularityOptimizer seurust/Seurat parity")

node1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
           1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 8, 8, 8, 9, 13,
           14, 14, 15, 15, 18, 18, 19, 20, 20, 22, 22, 23, 23, 23, 23, 23, 24,
           24, 24, 25, 26, 26, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32)
node2 <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 17, 19, 21, 31, 2, 3, 7, 13,
           17, 19, 21, 30, 3, 7, 8, 9, 13, 27, 28, 32, 7, 12, 13, 6, 10, 6, 10,
           16, 16, 30, 32, 33, 33, 33, 32, 33, 32, 33, 32, 33, 33, 32, 33, 32,
           33, 25, 27, 29, 32, 33, 25, 27, 31, 31, 29, 33, 33, 31, 33, 32, 33,
           32, 33, 32, 33, 33)
connections <- sparseMatrix(i = node2 + 1, j = node1 + 1, x = 1.0)

run_both <- function(...) {
  skip_if_no_seurust()
  cpp <- Seurat:::RunModularityClusteringCpp(SNN = connections, ...)
  rust <- seurust::RunModularityClusteringCpp(SNN = connections, ...)
  list(cpp = cpp, rust = rust)
}

test_that("Algorithm 1 parity", {
  out <- run_both(
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 1,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 564,
    printOutput = 0,
    edgefilename = ""
  )
  expect_equal(out$rust, out$cpp)
})

test_that("Algorithm 2 parity", {
  out <- run_both(
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 2,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 2,
    printOutput = 0,
    edgefilename = ""
  )
  expect_equal(out$rust, out$cpp)
})

test_that("Algorithm 3 parity", {
  out <- run_both(
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 3,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 42,
    printOutput = 0,
    edgefilename = ""
  )
  expect_equal(out$rust, out$cpp)
})

test_that("Algorithm 4 parity", {
  out <- run_both(
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 4,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 99,
    printOutput = 0,
    edgefilename = ""
  )
  expect_equal(out$rust, out$cpp)
})

test_that("Modularity function 2 parity", {
  out <- run_both(
    modularityFunction = 2,
    resolution = 0.5,
    algorithm = 1,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 7,
    printOutput = 0,
    edgefilename = ""
  )
  expect_equal(out$rust, out$cpp)
})

test_that("Multiple random starts parity", {
  out <- run_both(
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 1,
    nRandomStarts = 3,
    nIterations = 2,
    randomSeed = 123,
    printOutput = 0,
    edgefilename = ""
  )
  expect_equal(out$rust, out$cpp)
})

test_that("Edge file input parity", {
  skip_if_no_seurust()
  edgefile <- tempfile(fileext = ".txt")
  on.exit(unlink(edgefile), add = TRUE)
  nn <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  snn <- Seurat:::ComputeSNN(nn_ranked = nn, prune = 0)
  Seurat:::WriteEdgeFile(snn = snn, filename = edgefile, display_progress = FALSE)
  cpp <- Seurat:::RunModularityClusteringCpp(
    SNN = connections,
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 1,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 1,
    printOutput = FALSE,
    edgefilename = edgefile
  )
  rust <- seurust::RunModularityClusteringCpp(
    SNN = connections,
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 1,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 1,
    printOutput = FALSE,
    edgefilename = edgefile
  )
  expect_equal(rust, cpp)
})
