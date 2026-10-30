seurat_native <- function(name) {
  skip_if_not_installed("Seurat")
  getFromNamespace(x = name, ns = "Seurat")
}

seurat_has <- function(name) {
  skip_if_not_installed("Seurat")
  if (!exists(name, envir = asNamespace("Seurat"), inherits = FALSE)) {
    skip(paste("Installed Seurat does not provide", name))
  }
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

align_pc_signs <- function(a, b) {
  for (j in seq_len(ncol(a))) {
    if (sum(a[, j] * b[, j]) < 0) {
      a[, j] <- -a[, j]
    }
  }
  a
}
