seurat_native <- function(name) {
  skip_if_not_installed("Seurat")
  getFromNamespace(x = name, ns = "Seurat")
}

knn_graph <- function(n = 60, k = 6, seed = 1) {
  set.seed(seed)
  emb <- matrix(rnorm(n * 4), nrow = n)
  d <- as.matrix(dist(emb))
  idx <- t(apply(d, 1, function(r) order(r)[2:(k + 1)]))
  vals <- t(vapply(seq_len(n), function(i) d[i, idx[i, ]], numeric(k)))
  list(
    graph = Matrix::sparseMatrix(
      i = rep(seq_len(n), each = k),
      j = as.vector(t(idx)),
      x = as.vector(t(vals)),
      dims = c(n, n)
    ),
    idx = idx
  )
}

test_that("fast_dist accepts integer and double neighbor ids", {
  set.seed(3)
  x <- matrix(rnorm(40), nrow = 10)
  y <- matrix(rnorm(60), nrow = 15)
  n_int <- lapply(1:10, function(i) sample.int(15, 5))
  names(n_int) <- paste0("cell", 1:10)
  n_dbl <- lapply(n_int, as.double)

  got_int <- fast_dist(x, y, n_int)
  got_dbl <- fast_dist(x, y, n_dbl)
  expect_identical(got_int, got_dbl)
  expect_named(got_int, names(n_int))
  expect_equal(got_int[[2]], sqrt(colSums((t(y[n_int[[2]], ]) - x[2, ])^2)))

  expect_identical(fast_dist(x, y, n_int[1:3]), list())
  expect_error(fast_dist(x, y, c(list(16L), n_int[-1])), "outside")
  expect_error(fast_dist(x, y, c(list(NA_integer_), n_int[-1])), "outside")
})

test_that("fast_dist matches Seurat for integer and double ids", {
  seurat_fast_dist <- seurat_native("fast_dist")
  set.seed(4)
  x <- matrix(rnorm(80), nrow = 20)
  y <- matrix(rnorm(100), nrow = 25)
  n_int <- lapply(1:20, function(i) sample.int(25, 6))
  expect_equal(fast_dist(x, y, n_int), seurat_fast_dist(x, y, n_int), tolerance = 1e-12)
  n_dbl <- lapply(n_int, as.double)
  expect_equal(fast_dist(x, y, n_dbl), seurat_fast_dist(x, y, n_dbl), tolerance = 1e-12)
})

test_that("GraphToNeighborHelper reads neighbors from rows", {
  g <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3, 3),
    j = c(2, 3, 1, 3, 1, 2),
    x = c(0.5, 0.2, 0.9, 0.1, 0.3, 0.4),
    dims = c(3, 3)
  )
  out <- GraphToNeighborHelper(g)
  expect_equal(out[[1]], matrix(c(3, 3, 1, 2, 1, 2), nrow = 3))
  expect_equal(out[[2]], matrix(c(0.2, 0.1, 0.3, 0.5, 0.9, 0.4), nrow = 3))

  ragged <- Matrix::sparseMatrix(i = c(1, 1, 2), j = c(2, 3, 1), x = 1, dims = c(3, 3))
  expect_error(GraphToNeighborHelper(ragged), "equal number of neighbors")
})

test_that("GraphToNeighborHelper matches Seurat on an asymmetric kNN graph", {
  seurat_helper <- seurat_native("GraphToNeighborHelper")
  g <- knn_graph()$graph
  expect_false(Matrix::isSymmetric(g))
  expect_identical(GraphToNeighborHelper(g), seurat_helper(g))
})

test_that("RowMergeMatrices accepts both sparse layouts", {
  set.seed(5)
  m1 <- Matrix::rsparsematrix(8, 4, 0.4)
  m2 <- Matrix::rsparsematrix(6, 3, 0.4)
  rn1 <- paste0("g", 1:8)
  rn2 <- paste0("g", 5:10)
  all_rn <- union(rn1, rn2)
  from_csc <- RowMergeMatrices(m1, m2, rn1, rn2, all_rn)
  from_csr <- RowMergeMatrices(
    as(m1, "RsparseMatrix"), as(m2, "RsparseMatrix"), rn1, rn2, all_rn
  )
  expect_s4_class(from_csc, "dgCMatrix")
  expect_identical(from_csc, from_csr)
  expect_equal(dim(from_csc), c(10L, 7L))
  expect_equal(as.matrix(from_csc[5:8, 1:4]), as.matrix(m1[5:8, ]))
  expect_equal(as.matrix(from_csc[5:10, 5:7]), as.matrix(m2))

  seurat_merge <- seurat_native("RowMergeMatrices")
  expect_identical(
    from_csc,
    seurat_merge(as(m1, "RsparseMatrix"), as(m2, "RsparseMatrix"), rn1, rn2, all_rn)
  )
})

test_that("DirectSNNToFile returns a dgCMatrix and writes the edge file", {
  nn <- cbind(1:60, knn_graph()$idx)
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path))
  snn <- DirectSNNToFile(nn, prune = 1 / 15, display_progress = FALSE, filename = path)
  expect_s4_class(snn, "dgCMatrix")
  expect_identical(snn, ComputeSNN(nn, prune = 1 / 15))

  edges <- read.delim(path, header = FALSE)
  lower <- Matrix::summary(Matrix::tril(snn, k = -1))
  expect_equal(nrow(edges), nrow(lower))

  storage.mode(nn) <- "integer"
  expect_identical(ComputeSNN(nn, prune = 1 / 15), snn)
})

test_that("DirectSNNToFile and WriteEdgeFile match Seurat byte-for-byte", {
  seurat_direct <- seurat_native("DirectSNNToFile")
  nn <- cbind(1:60, knn_graph(seed = 7)$idx)
  f_seurat <- tempfile()
  f_rust <- tempfile()
  on.exit(unlink(c(f_seurat, f_rust)))
  expected <- seurat_direct(nn, 1 / 15, FALSE, f_seurat)
  got <- DirectSNNToFile(nn, 1 / 15, FALSE, f_rust)
  expect_equal(got, expected, tolerance = 1e-12)
  expect_identical(readLines(f_rust), readLines(f_seurat))

  WriteEdgeFile(expected, f_rust, FALSE)
  expect_identical(readLines(f_rust), readLines(f_seurat))
})

test_that("dense kernels match base R", {
  set.seed(6)
  a <- matrix(rnorm(120 * 70), nrow = 120)
  b <- matrix(rnorm(120 * 40), nrow = 120)
  expect_equal(FastCov(a), cov(a), tolerance = 1e-12)
  expect_true(isSymmetric(FastCov(a)))
  expect_equal(FastCov(a, center = FALSE), crossprod(a) / 119, tolerance = 1e-12)
  expect_equal(FastCovMats(a, b), cov(a, b), tolerance = 1e-12)
  expect_identical(FastRBind(a, a[1:5, ]), rbind(a, a[1:5, ]))
  expect_error(FastRBind(a, b), "columns")
  expect_error(FastCovMats(a, b[1:10, ]), "rows")

  ai <- matrix(sample.int(9, 60, replace = TRUE), nrow = 12)
  expect_equal(FastCov(ai), cov(ai), tolerance = 1e-12)
  expect_equal(RowVar(ai), apply(ai, 1, var), tolerance = 1e-12)
})

test_that("dense kernels match Seurat", {
  seurat_cov <- seurat_native("FastCov")
  seurat_cov_mats <- seurat_native("FastCovMats")
  seurat_rbind <- seurat_native("FastRBind")
  set.seed(8)
  a <- matrix(rnorm(300 * 90), nrow = 300)
  b <- matrix(rnorm(300 * 50), nrow = 300)
  expect_equal(FastCov(a), seurat_cov(a), tolerance = 1e-12)
  expect_equal(FastCovMats(a, b), seurat_cov_mats(a, b), tolerance = 1e-12)
  expect_identical(FastRBind(a, a), seurat_rbind(a, a))
})

test_that("numeric vector arguments accept integers", {
  set.seed(9)
  mat <- Matrix::rsparsematrix(20, 30, 0.3)
  mu <- sample.int(3, 20, replace = TRUE)
  expect_equal(SparseRowVar2(mat, mu, FALSE), SparseRowVar2(mat, as.double(mu), FALSE))
})
