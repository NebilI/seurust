system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")

suppressPackageStartupMessages({
  devtools::load_all(recompile = FALSE, quiet = TRUE)
  library(seurust)
  library(Matrix)
})

cat("==> Parity: sparse row stats...\n")
m <- sparseMatrix(
  i = c(1L, 3L, 2L, 3L),
  p = c(0L, 2L, 3L, 4L),
  x = c(1, 2, 3, 4),
  dims = c(3L, 3L)
)
x <- slot(m, "x")
i <- slot(m, "i")
stopifnot(all.equal(
  Seurat:::row_sum_dgcmatrix(x, i, nrow(m), ncol(m)),
  seurust::row_sum_dgcmatrix(x, i, nrow(m), ncol(m))
))
cat("Row stats OK\n")

cat("==> Parity: log normalization...\n")
mat <- as(matrix(1:16, ncol = 4, nrow = 4), "sparseMatrix")
cpp <- Seurat:::LogNorm(mat, 1e4, display_progress = FALSE)
rust <- seurust::LogNorm(mat, 1e4, display_progress = FALSE)
stopifnot(all.equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10))
cat("LogNorm OK\n")

cat("==> Parity: dense covariance...\n")
set.seed(42)
mat <- replicate(10, rchisq(10, 4))
stopifnot(all.equal(Seurat:::FastCov(mat), seurust::FastCov(mat)))
cat("FastCov OK\n")

cat("==> Parity: ComputeSNN...\n")
nn <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
cpp <- Seurat:::ComputeSNN(nn, 0.01)
rust <- seurust::ComputeSNN(nn, 0.01)
stopifnot(all.equal(as.matrix(cpp), as.matrix(rust), tolerance = 1e-10))
cat("ComputeSNN OK\n")

run_tests <- function(path) {
  res <- testthat::test_file(path, reporter = "summary", stop_on_failure = TRUE)
  invisible(res)
}

cat("==> Running modularity optimizer tests...\n")
run_tests("tests/testthat/test_modularity_optimizer.R")
run_tests("tests/testthat/test_rust_cpp_parity_modularity.R")

cat("==> Timing: Seurat vs seurust (ratio > 1 => Rust faster)...\n")
run_tests("tests/testthat/test_rust_cpp_timing_modularity.R")
run_tests("tests/testthat/test_rust_cpp_timing.R")

cat("All parity checks passed.\n")
