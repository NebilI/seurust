test_that("LogNorm matches Seurat 5.6 slot API", {
  seurat_log <- seurat_has("LogNorm")
  set.seed(11)
  mat <- as(Matrix::rsparsematrix(40, 30, 0.2, rand.x = function(n) rpois(n, 4) + 1), "dgCMatrix")
  rust <- seurust::LogNorm(mat, scale_factor = 1e4, display_progress = FALSE)
  cpp_x <- seurat_log(mat@x, mat@p, 1e4, 1L, FALSE)
  expect_s4_class(rust, "dgCMatrix")
  expect_equal(as.numeric(rust@x), as.numeric(cpp_x), tolerance = 1e-12)
})

test_that("ComputeSNN matches Seurat with nthreads", {
  seurat_snn <- seurat_has("ComputeSNN")
  nn <- cbind(1:40, knn_graph(n = 40, k = 5, seed = 12)$idx)
  rust <- seurust::ComputeSNN(nn, prune = 1 / 15, nthreads = 1L)
  cpp <- seurat_snn(nn, 1 / 15, 1L)
  expect_equal(rust, cpp, tolerance = 1e-12)
})

test_that("FastDenseRowScale matches Seurat", {
  seurat_scale <- seurat_has("FastDenseRowScale")
  set.seed(13)
  mat <- matrix(rnorm(80 * 25), nrow = 80)
  expect_equal(
    seurust::FastDenseRowScale(mat, display_progress = FALSE),
    seurat_scale(mat, display_progress = FALSE),
    tolerance = 1e-10
  )
  feats <- as.integer(c(0L, 3L, 10L, 40L))
  expect_equal(
    seurust::FastDenseRowScale(mat, features = feats, display_progress = FALSE),
    seurat_scale(mat, features = feats, display_progress = FALSE),
    tolerance = 1e-10
  )
})

test_that("SparseRowMeanVar and legacy match Seurat", {
  seurat_mv <- seurat_has("SparseRowMeanVar")
  seurat_legacy <- seurat_has("SparseRowMeanVarLegacy")
  set.seed(14)
  mat <- as(Matrix::rsparsematrix(50, 40, 0.15, rand.x = function(n) rpois(n, 3) + 1), "dgCMatrix")
  rust <- seurust::SparseRowMeanVar(mat@x, mat@i, mat@p, nrow(mat), ncol(mat), 1L, FALSE)
  cpp <- seurat_mv(mat@x, mat@i, mat@p, nrow(mat), ncol(mat), 1L, FALSE)
  expect_equal(as.numeric(rust$mean), as.numeric(cpp$mean), tolerance = 1e-12)
  expect_equal(as.numeric(rust$variance), as.numeric(cpp$variance), tolerance = 1e-12)

  rust_l <- seurust::SparseRowMeanVarLegacy(mat@x, mat@i, mat@p, nrow(mat), ncol(mat), FALSE)
  cpp_l <- seurat_legacy(mat@x, mat@i, mat@p, nrow(mat), ncol(mat), FALSE)
  expect_equal(as.numeric(rust_l$mean), as.numeric(cpp_l$mean), tolerance = 1e-12)
  expect_equal(as.numeric(rust_l$variance), as.numeric(cpp_l$variance), tolerance = 1e-12)
})

test_that("SparseRowVarStd matches Seurat 5.6 slot API", {
  seurat_std <- seurat_has("SparseRowVarStd")
  set.seed(15)
  mat <- as(Matrix::rsparsematrix(30, 50, 0.2, rand.x = function(n) rpois(n, 2) + 1), "dgCMatrix")
  mu <- as.numeric(Matrix::rowMeans(mat))
  sd <- sqrt(pmax(mu, 1e-8))
  rust <- seurust::SparseRowVarStd(mat, mu, sd, 10, FALSE)
  cpp <- seurat_std(mat@x, mat@i, mat@p, mu, sd, nrow(mat), ncol(mat), 10, 1L, FALSE)
  expect_equal(as.numeric(rust), as.numeric(cpp), tolerance = 1e-12)
})

test_that("SparseRowVarStdLegacyRows matches Seurat", {
  seurat_rows <- seurat_has("SparseRowVarStdLegacyRows")
  set.seed(16)
  mat <- as(Matrix::rsparsematrix(40, 35, 0.2, rand.x = function(n) rpois(n, 2) + 1), "dgCMatrix")
  mu <- as.numeric(Matrix::rowMeans(mat))
  sd <- sqrt(pmax(mu, 1e-8))
  rows_use <- as.integer(c(0L, 2L, 7L, 19L, 33L))
  rust <- seurust::SparseRowVarStdLegacyRows(
    mat@x, mat@i, mat@p, mu, sd, 10, nrow(mat), ncol(mat), rows_use, 1L, FALSE
  )
  cpp <- seurat_rows(
    mat@x, mat@i, mat@p, mu, sd, 10, nrow(mat), ncol(mat), rows_use, 1L, FALSE
  )
  expect_equal(as.numeric(rust), as.numeric(cpp), tolerance = 1e-12)
})

test_that("FindAllMarkersSparseFoldChangeStats matches Seurat", {
  seurat_fc <- seurat_has("FindAllMarkersSparseFoldChangeStats")
  set.seed(17)
  mat <- as(Matrix::rsparsematrix(25, 60, 0.2, rand.x = function(n) rexp(n, 1)), "dgCMatrix")
  groups <- as.integer(rep(1:3, length.out = ncol(mat)))
  rust <- seurust::FindAllMarkersSparseFoldChangeStats(
    mat@x, mat@i, mat@p, nrow(mat), ncol(mat), groups, 3L, TRUE, 1L
  )
  cpp <- seurat_fc(mat@x, mat@i, mat@p, nrow(mat), ncol(mat), groups, 3L, TRUE, 1L)
  expect_equal(as.numeric(rust$group_sum), as.numeric(cpp$group_sum), tolerance = 1e-10)
  expect_equal(as.numeric(rust$rest_sum), as.numeric(cpp$rest_sum), tolerance = 1e-10)
  expect_equal(as.integer(rust$detected), as.integer(cpp$detected))
  expect_equal(as.integer(rust$total_detected), as.integer(cpp$total_detected))
})

test_that("CcaCrossprodMultiply matches Seurat", {
  seurat_cca <- seurat_has("CcaCrossprodMultiply")
  set.seed(18)
  left <- matrix(rnorm(40 * 12), nrow = 40)
  right <- matrix(rnorm(40 * 9), nrow = 40)
  x <- rnorm(9)
  expect_equal(
    as.numeric(seurust::CcaCrossprodMultiply(left, right, x)),
    as.numeric(seurat_cca(left, right, x)),
    tolerance = 1e-10
  )
})

test_that("CountAnchorSharedNeighbors matches Seurat", {
  seurat_count <- seurat_has("CountAnchorSharedNeighbors")
  set.seed(19)
  n_a <- 20L
  n_b <- 18L
  k <- 5L
  mk <- function(nrow, ncol, nmax) {
    m <- matrix(sample.int(nmax, nrow * ncol, replace = TRUE), nrow = nrow)
    storage.mode(m) <- "double"
    m
  }
  aa <- mk(n_a, k, n_a)
  ab <- mk(n_a, k, n_b)
  ba <- mk(n_b, k, n_a)
  bb <- mk(n_b, k, n_b)
  anchors <- 8L
  cell1 <- sample.int(n_a, anchors, replace = TRUE)
  cell2 <- sample.int(n_b, anchors, replace = TRUE)
  rust <- seurust::CountAnchorSharedNeighbors(aa, ab, ba, bb, cell1, cell2, n_a, k, 1L)
  cpp <- seurat_count(aa, ab, ba, bb, cell1, cell2, n_a, k, 1L)
  expect_equal(as.integer(rust), as.integer(cpp))
})

test_that("EigenGramPCA matches Seurat up to loading sign", {
  seurat_pca <- seurat_has("EigenGramPCA")
  set.seed(20)
  object <- matrix(rnorm(30 * 80), nrow = 30)
  rust <- seurust::EigenGramPCA(object, npcs = 5L, weight_by_var = TRUE, nthreads = 1L)
  cpp <- seurat_pca(object, 5L, TRUE, 1L)
  rust$loadings <- align_pc_signs(rust$loadings, cpp$loadings)
  rust$embeddings <- align_pc_signs(rust$embeddings, cpp$embeddings)
  expect_equal(rust$loadings, cpp$loadings, tolerance = 1e-8)
  expect_equal(rust$embeddings, cpp$embeddings, tolerance = 1e-8)
  expect_equal(as.numeric(rust$sdev), as.numeric(cpp$sdev), tolerance = 1e-8)
})

test_that("SCT residual kernels match Seurat", {
  seurat_stats <- seurat_has("SCTResidualStatsAndCorrected")
  seurat_mat <- seurat_has("SCTPearsonResidualMatrix")
  set.seed(21)
  mat <- as(Matrix::rsparsematrix(20, 40, 0.25, rand.x = function(n) rpois(n, 3)), "dgCMatrix")
  theta <- runif(nrow(mat), 1, 8)
  intercept <- rnorm(nrow(mat), -1, 0.2)
  slope <- rep(0.9, nrow(mat))
  log_umi <- log1p(Matrix::colSums(mat))
  clip <-  sqrt(ncol(mat))
  rust <- seurust::SCTResidualStatsAndCorrected(
    mat@x, mat@i, mat@p, nrow(mat), ncol(mat),
    theta, intercept, slope, log_umi,
    median(log_umi), 0.04, -clip, clip, 1L, TRUE
  )
  cpp <- seurat_stats(
    mat@x, mat@i, mat@p, nrow(mat), ncol(mat),
    theta, intercept, slope, log_umi,
    median(log_umi), 0.04, -clip, clip, 1L, TRUE
  )
  expect_equal(as.numeric(rust$residual_mean), as.numeric(cpp$residual_mean), tolerance = 1e-10)
  expect_equal(as.numeric(rust$residual_variance), as.numeric(cpp$residual_variance), tolerance = 1e-10)
  expect_equal(as.matrix(rust$corrected), as.matrix(cpp$corrected), tolerance = 1e-10)

  feats <- as.integer(c(0L, 3L, 8L, 15L))
  rust_y <- seurust::SCTPearsonResidualMatrix(
    mat@x, mat@i, mat@p, nrow(mat), ncol(mat),
    theta, intercept, slope, log_umi, feats, 0.04, -clip, clip, TRUE, 1L, FALSE
  )
  cpp_y <- seurat_mat(
    mat@x, mat@i, mat@p, nrow(mat), ncol(mat),
    theta, intercept, slope, log_umi, feats, 0.04, -clip, clip, TRUE, 1L, FALSE
  )
  expect_equal(rust_y, cpp_y, tolerance = 1e-10)
})

test_that("RunModularityClusteringCpp_multi matches repeated single runs", {
  skip_if_not_installed("Seurat")
  nn <- cbind(1:40, knn_graph(n = 40, k = 6, seed = 22)$idx)
  snn <- seurust::ComputeSNN(nn, prune = 0.01)
  resolutions <- c(0.4, 0.8, 1.2)
  multi <- seurust::RunModularityClusteringCpp_multi(
    snn, 1L, resolutions, 1L, 1L, 5L, 7L, FALSE, "", 1L
  )
  expect_length(multi, 3L)
  for (i in seq_along(resolutions)) {
    single <- seurust::RunModularityClusteringCpp(
      snn, 1L, resolutions[[i]], 1L, 1L, 5L, 7L, FALSE, "", 1L
    )
    expect_identical(as.integer(multi[[i]]), as.integer(single))
  }
})
