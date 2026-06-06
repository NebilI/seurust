# Timing comparison: Seurat (C++) vs SeuratRust modularity clustering.
context("ModularityOptimizer SeuratRust/Seurat timing")

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

modularity_args <- list(
  modularityFunction = 1L,
  resolution = 1.0,
  algorithm = 3L,
  nRandomStarts = 5L,
  nIterations = 50L,
  randomSeed = 42L,
  printOutput = FALSE,
  edgefilename = ""
)

test_that("Modularity clustering timing", {
  skip_if_no_seuratrust()
  run_cpp <- function() {
    do.call(Seurat:::RunModularityClusteringCpp, c(list(SNN = connections), modularity_args))
  }
  run_rust <- function() {
    do.call(SeuratRust::RunModularityClusteringCpp, c(list(SNN = connections), modularity_args))
  }
  out_cpp <- run_cpp()
  out_rust <- run_rust()
  expect_equal(out_rust, out_cpp)
  bench <- benchmark_rust_cpp(
    cpp_fn = run_cpp,
    rust_fn = run_rust,
    n_warmup = 2L,
    n_reps = 10L
  )
  expect_timing_report(bench, "Modularity clustering")
  if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
    expect_rust_faster(bench, "Modularity clustering", tolerance = 0.95)
  }
})
