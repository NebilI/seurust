# High-level R wrappers matching Seurat's RcppExports API (Rust backend).

#' Coerce to a double-precision base matrix
#'
#' Seurat's Rcpp signatures (`NumericMatrix`, `Eigen::MatrixXd`) silently
#' accept integer and logical matrices; the extendr entry points only accept
#' doubles, so coerce before crossing into Rust.
#' @keywords internal
#' @noRd
AsDoubleMatrix <- function(x) {
  if (!is.matrix(x = x)) {
    x <- as.matrix(x = x)
  }
  if (!is.double(x = x)) {
    storage.mode(x = x) <- "double"
  }
  x
}

#' Coerce any matrix to a double-precision general sparse matrix
#' @keywords internal
#' @noRd
AsGeneralSparse <- function(mat, layout) {
  as(object = as(object = as(object = mat, Class = "dMatrix"), Class = "generalMatrix"), Class = layout)
}

#' Extract dgCMatrix slots for extendr calls
#' @keywords internal
#' @noRd
CscSlots <- function(mat) {
  if (!inherits(x = mat, what = "dgCMatrix")) {
    mat <- AsGeneralSparse(mat = mat, layout = "CsparseMatrix")
  }
  list(
    x = slot(object = mat, name = "x"),
    i = slot(object = mat, name = "i"),
    p = slot(object = mat, name = "p"),
    nrows = nrow(x = mat),
    ncols = ncol(x = mat)
  )
}

#' Extract dgRMatrix slots for extendr calls
#' @keywords internal
#' @noRd
CsrSlots <- function(mat) {
  if (!inherits(x = mat, what = "dgRMatrix")) {
    mat <- AsGeneralSparse(mat = mat, layout = "RsparseMatrix")
  }
  list(
    x = slot(object = mat, name = "x"),
    j = slot(object = mat, name = "j"),
    p = slot(object = mat, name = "p"),
    nrows = nrow(x = mat),
    ncols = ncol(x = mat)
  )
}

#' Reconstruct a dgCMatrix from an extendr slot list
#' @keywords internal
#' @noRd
CscFromList <- function(slots) {
  methods::new(
    Class = "dgCMatrix",
    x = slots$x,
    i = slots$i,
    p = slots$p,
    Dim = as.integer(slots$Dim)
  )
}

#' Log-normalize a sparse count matrix
#'
#' Rust reimplementation of Seurat's `LogNorm`. Divides each column by its
#' sum, multiplies by `scale_factor`, then applies `log1p`.
#'
#' @param data A `dgCMatrix` (or coercible) of counts (features x cells).
#' @param scale_factor Numeric scale factor (Seurat default is `1e4`).
#' @param display_progress Logical; show a progress indicator when supported.
#'
#' @return A log-normalized `dgCMatrix` with the same dimensions as `data`.
#'
#' @examples
#' library(Matrix)
#' mat <- sparseMatrix(i = c(1, 3, 2), j = c(1, 2, 3), x = 1:3, dims = c(3, 3))
#' out <- LogNorm(mat, scale_factor = 1e4, display_progress = FALSE)
#' out[1, 1]
#'
#' @export
LogNorm <- function(data, scale_factor, display_progress = TRUE) {
  if (!inherits(x = data, what = "dgCMatrix")) {
    data <- AsGeneralSparse(mat = data, layout = "CsparseMatrix")
  }
  s <- CscSlots(mat = data)
  slot(object = data, name = "x") <- log_norm(
    x = s$x,
    i = s$i,
    p = s$p,
    nrows = s$nrows,
    ncols = s$ncols,
    scale_factor = scale_factor,
    display_progress = display_progress
  )
  data
}

#' Standardize columns of a dense matrix
#'
#' @param mat Numeric matrix.
#' @param display_progress Logical; show progress when supported.
#'
#' @return A numeric matrix with standardized columns.
#' @export
Standardize <- function(mat, display_progress = TRUE, nthreads = 1L) {
  standardize(
    mat = AsDoubleMatrix(x = mat),
    display_progress = display_progress,
    `_nthreads` = as.integer(nthreads)
  )
}

#' Fast covariance of a dense matrix
#'
#' @param mat Numeric matrix.
#' @param center Logical; center columns before computing covariance.
#'
#' @return A symmetric numeric covariance matrix.
#' @export
FastCov <- function(mat, center = TRUE) {
  fast_cov(mat = AsDoubleMatrix(x = mat), center = center)
}

#' Fast cross-covariance of two dense matrices
#'
#' @param mat1,mat2 Numeric matrices with the same number of rows.
#' @param center Logical; center columns before computing covariance.
#'
#' @return A numeric cross-covariance matrix.
#' @export
FastCovMats <- function(mat1, mat2, center = TRUE) {
  fast_cov_mats(
    mat1 = AsDoubleMatrix(x = mat1),
    mat2 = AsDoubleMatrix(x = mat2),
    center = center
  )
}

#' Fast row-bind of two dense matrices
#'
#' @param mat1,mat2 Numeric matrices with the same number of columns.
#'
#' @return A numeric matrix stacking `mat1` and `mat2` by row.
#' @export
FastRBind <- function(mat1, mat2) {
  fast_rbind(mat1 = AsDoubleMatrix(x = mat1), mat2 = AsDoubleMatrix(x = mat2))
}

#' Row variances of a dense matrix
#'
#' @param x Numeric matrix.
#'
#' @return Numeric vector of row variances.
#' @export
RowVar <- function(x) {
  row_var(mat = AsDoubleMatrix(x = x))
}

#' Merge two sparse matrices by shared row names
#'
#' Seurat's C++ entry point only accepts row-compressed (`RsparseMatrix`)
#' inputs and errors on a `dgCMatrix`. This version accepts either layout (or
#' any coercible matrix) and converts to row-compressed form first; for
#' `RsparseMatrix` inputs the result is identical to Seurat's.
#'
#' @param mat1,mat2 Sparse matrices, preferably row-compressed
#'   (`dgRMatrix`); other matrix classes are coerced.
#' @param mat1_rownames,mat2_rownames Character vectors of row names.
#' @param all_rownames Character vector of the union of row names.
#'
#' @return A `dgCMatrix` containing the merged rows.
#' @export
RowMergeMatrices <- function(mat1, mat2, mat1_rownames, mat2_rownames, all_rownames) {
  s1 <- CsrSlots(mat = mat1)
  s2 <- CsrSlots(mat = mat2)
  row_merge_matrices(
    x1 = s1$x, j1 = s1$j, p1 = s1$p, nrows1 = s1$nrows, ncols1 = s1$ncols,
    x2 = s2$x, j2 = s2$j, p2 = s2$p, nrows2 = s2$nrows, ncols2 = s2$ncols,
    mat1_rownames = mat1_rownames,
    mat2_rownames = mat2_rownames,
    all_rownames = all_rownames
  )
}

#' Replace columns of a sparse matrix
#'
#' @param mat A `dgCMatrix`.
#' @param col_idx Integer indices of columns to replace (0- or 1-based as in Seurat).
#' @param replacement A `dgCMatrix` whose columns replace those in `mat`.
#'
#' @return A `dgCMatrix` with the specified columns replaced.
#' @export
ReplaceColsC <- function(mat, col_idx, replacement) {
  s <- CscSlots(mat = mat)
  r <- CscSlots(mat = replacement)
  replace_cols(
    x = s$x, i = s$i, p = s$p, nrows = s$nrows, ncols = s$ncols,
    col_idx = as.double(x = col_idx),
    rx = r$x, ri = r$i, rp = r$p, rnrows = r$nrows, rncols = r$ncols
  )
}

#' Convert a sparse graph to neighbor index lists
#'
#' Reads each row of `mat` (cell `k`'s neighbors are the non-zero columns of
#' row `k`), orders them by increasing value, and returns them as matrices,
#' exactly like Seurat's `GraphToNeighborHelper`. Every row must have the same
#' number of non-zero entries.
#'
#' @param mat A sparse adjacency/`Graph` matrix (`dgCMatrix` or coercible),
#'   cells x cells.
#'
#' @return An unnamed list of two `nrow(mat)` x `k` numeric matrices: 1-based
#'   neighbor indices and the corresponding distances, both sorted by distance.
#' @export
GraphToNeighborHelper <- function(mat) {
  s <- CscSlots(mat = mat)
  graph_to_neighbor_helper(
    x = s$x, i = s$i, p = s$p, nrows = s$nrows, ncols = s$ncols
  )
}

#' Fast exp-mean of sparse matrix rows
#'
#' @param mat A `dgCMatrix`.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Numeric vector of row-wise exp-means.
#' @export
FastExpMean <- function(mat, display_progress) {
  s <- CscSlots(mat = mat)
  fast_exp_mean(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

#' Sparse row variances
#'
#' @param mat A `dgCMatrix`.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Numeric vector of row variances.
#' @export
SparseRowVar <- function(mat, display_progress) {
  s <- CscSlots(mat = mat)
  sparse_row_var(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

#' Sparse row variances given means
#'
#' @param mat A `dgCMatrix`.
#' @param mu Numeric vector of row means.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Numeric vector of row variances.
#' @export
SparseRowVar2 <- function(mat, mu, display_progress) {
  s <- CscSlots(mat = mat)
  sparse_row_var2(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    mu = as.double(x = mu),
    display_progress = display_progress
  )
}

#' Sparse row variances with clipping
#'
#' @param mat A `dgCMatrix`.
#' @param mu Numeric vector of row means.
#' @param sd Numeric vector of row standard deviations.
#' @param vmax Maximum standardized value.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Numeric vector of clipped row variances.
#' @export
SparseRowVarStd <- function(mat, mu, sd, vmax, display_progress) {
  s <- CscSlots(mat = mat)
  sparse_row_var_std(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    mu = as.double(x = mu), sd = as.double(x = sd), vmax = vmax,
    display_progress = display_progress
  )
}

#' Fast log variance-to-mean ratio
#'
#' @param mat A `dgCMatrix`.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Numeric vector of log VMR values.
#' @export
FastLogVMR <- function(mat, display_progress) {
  s <- CscSlots(mat = mat)
  fast_log_vmr(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

#' Scale and/or center sparse matrix rows
#'
#' @param mat A `dgCMatrix`.
#' @param scale Logical; scale rows.
#' @param center Logical; center rows.
#' @param scale_max Maximum absolute scaled value.
#' @param display_progress Logical; show progress when supported.
#'
#' @return A dense numeric matrix of scaled values.
#' @export
FastSparseRowScale <- function(
    mat,
    scale = TRUE,
    center = TRUE,
    scale_max = 10,
    display_progress = TRUE,
    features = integer(),
    nthreads = 1L) {
  s <- CscSlots(mat = mat)
  fast_sparse_row_scale(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    features = as.integer(x = features),
    scale = scale, center = center, scale_max = scale_max,
    display_progress = display_progress
  )
}

#' Scale sparse matrix rows with known statistics
#'
#' @inheritParams FastSparseRowScale
#' @param mu Numeric vector of row means.
#' @param sigma Numeric vector of row standard deviations.
#'
#' @return A dense numeric matrix of scaled values.
#' @export
FastSparseRowScaleWithKnownStats <- function(mat, mu, sigma, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  s <- CscSlots(mat = mat)
  fast_sparse_row_scale_with_known_stats(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    mu = as.double(x = mu), sigma = as.double(x = sigma),
    scale = scale, center = center, scale_max = scale_max,
    display_progress = display_progress
  )
}

#' Downsample UMIs globally
#'
#' @param data A `dgCMatrix` of counts.
#' @param sample_val Target sampling value.
#' @param upsample Logical; allow upsampling.
#' @param display_progress Logical; show progress when supported.
#'
#' @return A sampled `dgCMatrix`.
#' @export
RunUMISampling <- function(data, sample_val, upsample = FALSE, display_progress = TRUE) {
  s <- CscSlots(mat = data)
  out <- run_umi_sampling(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    sample_val = sample_val, upsample = upsample,
    `_display_progress` = display_progress
  )
  CscFromList(out)
}

#' Downsample UMIs per cell
#'
#' @inheritParams RunUMISampling
#'
#' @return A sampled `dgCMatrix`.
#' @export
RunUMISamplingPerCell <- function(data, sample_val, upsample = FALSE, display_progress = TRUE) {
  s <- CscSlots(mat = data)
  out <- run_umi_sampling_per_cell(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    sample_val = as.double(x = sample_val), upsample = upsample,
    `_display_progress` = display_progress
  )
  CscFromList(out)
}

#' Build a shared nearest-neighbor (SNN) graph
#'
#' @param nn_ranked Matrix of 1-based neighbor indices (cells x k); integer
#'   or double.
#' @param prune Numeric pruning threshold.
#'
#' @return A sparse SNN graph as a `dgCMatrix`.
#' @export
ComputeSNN <- function(nn_ranked, prune, nthreads = 1L) {
  compute_snn(
    nn_ranked = AsDoubleMatrix(x = nn_ranked),
    prune = prune,
    `_nthreads` = as.integer(nthreads)
  )
}

#' Integrate expression using anchor weights
#'
#' @param integration_matrix Sparse integration matrix.
#' @param weights Sparse weight matrix.
#' @param expression_cells2 Sparse expression matrix for dataset 2.
#'
#' @return A corrected `dgCMatrix`.
#' @export
IntegrateDataC <- function(integration_matrix, weights, expression_cells2) {
  im <- CscSlots(integration_matrix)
  w <- CscSlots(weights)
  ex <- CscSlots(expression_cells2)
  integrate_data(
    ix = im$x, ii = im$i, ip = im$p, inrows = im$nrows, incols = im$ncols,
    wx = w$x, wi = w$i, wp = w$p, wnrows = w$nrows, wncols = w$ncols,
    ex = ex$x, ei = ex$i, ep = ex$p, enrows = ex$nrows, encols = ex$ncols
  )
}

#' Find integration anchor weights
#'
#' @param cells2 Character vector of cell names in dataset 2.
#' @param distances Numeric distance matrix/vector as used by Seurat.
#' @param anchor_cells2 Character vector of anchor cells.
#' @param integration_matrix_rownames Character vector of integration matrix rows.
#' @param cell_index Integer cell indices.
#' @param anchor_score Numeric anchor scores.
#' @param min_dist,sd Numeric distance and bandwidth parameters.
#' @param display_progress Logical; show progress when supported.
#'
#' @return A sparse weight `dgCMatrix`.
#' @export
FindWeightsC <- function(cells2, distances, anchor_cells2, integration_matrix_rownames,
                         cell_index, anchor_score, min_dist, sd, display_progress,
                         nthreads = 1L) {
  find_weights(
    cells2 = as.double(x = cells2),
    distances = AsDoubleMatrix(x = distances),
    anchor_cells2 = anchor_cells2,
    integration_matrix_rownames = integration_matrix_rownames,
    cell_index = AsDoubleMatrix(x = cell_index),
    anchor_score = as.double(x = anchor_score),
    min_dist = min_dist,
    sd = sd,
    display_progress = display_progress,
    `_nthreads` = as.integer(nthreads)
  )
}

#' Score query cells using an SNN graph
#'
#' @param snn Sparse SNN graph.
#' @param query_pca Query PCA embeddings.
#' @param query_dists Query distances.
#' @param corrected_nns Corrected nearest neighbors.
#' @param k_snn Integer number of SNN neighbors.
#' @param subtract_first_nn Logical; subtract first nearest neighbor.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Numeric scores matching Seurat's helper.
#' @export
ScoreHelper <- function(snn, query_pca, query_dists, corrected_nns, k_snn,
                        subtract_first_nn, display_progress, nthreads = 1L) {
  s <- CscSlots(snn)
  score_helper(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    query_pca = AsDoubleMatrix(x = query_pca),
    query_dists = AsDoubleMatrix(x = query_dists),
    corrected_nns = AsDoubleMatrix(x = corrected_nns),
    k_snn = k_snn,
    subtract_first_nn = subtract_first_nn,
    display_progress = display_progress,
    `_nthreads` = as.integer(nthreads)
  )
}

#' Write SNN edges to a file
#'
#' @param snn Sparse SNN graph.
#' @param filename Output path.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Invisibly returns `NULL`.
#' @export
WriteEdgeFile <- function(snn, filename, display_progress) {
  s <- CscSlots(mat = snn)
  invisible(write_edge_file(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    filename = filename,
    display_progress = display_progress
  ))
}

#' Build an SNN graph and write it directly to a file
#'
#' Computes the SNN graph as [ComputeSNN()] does and writes its lower
#' triangle to `filename` in the same tab-separated format as
#' [WriteEdgeFile()].
#'
#' @param nn_ranked Matrix of 1-based neighbor indices (cells x k); integer
#'   or double.
#' @param prune Numeric pruning threshold.
#' @param display_progress Logical; show progress when supported.
#' @param filename Output path.
#'
#' @return The SNN graph as a `dgCMatrix`, as returned by Seurat.
#' @export
DirectSNNToFile <- function(nn_ranked, prune, display_progress, filename) {
  direct_snn_to_file(
    nn_ranked = AsDoubleMatrix(x = nn_ranked),
    prune = prune,
    display_progress = display_progress,
    filename = filename
  )
}

#' Smallest non-zero SNN distances
#'
#' @param snn Sparse SNN graph.
#' @param mat Numeric matrix of distances/embeddings as used by Seurat.
#' @param n Integer number of neighbors.
#' @param nearest_dist Numeric vector of nearest distances.
#'
#' @return Numeric vector of smallest non-zero distances.
#' @export
SNN_SmallestNonzero_Dist <- function(snn, mat, n, nearest_dist) {
  s <- CscSlots(snn)
  snn_smallest_nonzero_dist(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    mat = AsDoubleMatrix(x = mat), n = n, nearest_dist = as.double(x = nearest_dist)
  )
}

#' Run modularity clustering on an SNN graph
#'
#' Rust/C++ bridge matching Seurat's `RunModularityClusteringCpp`.
#'
#' @param SNN Sparse SNN graph (`dgCMatrix`).
#' @param modularityFunction Integer modularity function id.
#' @param resolution Numeric clustering resolution.
#' @param algorithm Integer algorithm id.
#' @param nRandomStarts Integer number of random starts.
#' @param nIterations Integer number of iterations.
#' @param randomSeed Integer random seed.
#' @param printOutput Logical; print progress.
#' @param edgefilename Optional edge list file path.
#'
#' @return Integer cluster assignments.
#' @export
RunModularityClusteringCpp <- function(
    SNN,
    modularityFunction = 1,
    resolution = 1.0,
    algorithm = 1,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 0,
    printOutput = FALSE,
    edgefilename = "",
    nThreads = 1L) {
  s <- CscSlots(mat = SNN)
  run_modularity_clustering(
    x = s$x,
    i = s$i,
    p = s$p,
    nrows = s$nrows,
    ncols = s$ncols,
    modularity_function = modularityFunction,
    resolution = resolution,
    algorithm = algorithm,
    n_random_starts = nRandomStarts,
    n_iterations = nIterations,
    random_seed = as.integer(randomSeed),
    print_output = as.logical(printOutput),
    edgefilename = edgefilename,
    n_threads = as.integer(nThreads)
  )
}

#' Run modularity clustering at multiple resolutions
#'
#' Rust reimplementation of Seurat 5.6 `RunModularityClusteringCpp_multi`.
#' Builds the network once and clusters at each resolution.
#'
#' @inheritParams RunModularityClusteringCpp
#' @param resolutions Numeric vector of clustering resolutions.
#'
#' @return A list of integer cluster assignment vectors, one per resolution.
#' @export
RunModularityClusteringCpp_multi <- function(
    SNN,
    modularityFunction = 1,
    resolutions,
    algorithm = 1,
    nRandomStarts = 1,
    nIterations = 1,
    randomSeed = 0,
    printOutput = FALSE,
    edgefilename = "",
    nThreads = 1L) {
  s <- CscSlots(mat = SNN)
  run_modularity_clustering_multi(
    x = s$x,
    i = s$i,
    p = s$p,
    nrows = s$nrows,
    ncols = s$ncols,
    modularity_function = modularityFunction,
    resolutions = as.double(x = resolutions),
    algorithm = algorithm,
    n_random_starts = nRandomStarts,
    n_iterations = nIterations,
    random_seed = as.integer(randomSeed),
    print_output = as.logical(printOutput),
    edgefilename = edgefilename,
    n_threads = as.integer(nThreads)
  )
}

#' Scale and/or center dense matrix rows
#'
#' Rust reimplementation of Seurat 5.6 `FastDenseRowScale`.
#'
#' @param mat Numeric matrix (features x cells).
#' @param features Optional 0-based integer row subset. Empty means all rows.
#' @param scale Logical; scale rows.
#' @param center Logical; center rows.
#' @param scale_max Maximum scaled value.
#' @param nthreads Integer thread hint (Rust uses the rayon pool).
#' @param display_progress Logical; show progress when supported.
#'
#' @return A dense numeric matrix of scaled values.
#' @export
FastDenseRowScale <- function(
    mat,
    features = integer(),
    scale = TRUE,
    center = TRUE,
    scale_max = 10,
    nthreads = 1L,
    display_progress = FALSE) {
  fast_dense_row_scale(
    mat = AsDoubleMatrix(x = mat),
    features = as.integer(x = features),
    scale = scale,
    center = center,
    scale_max = scale_max,
    `_nthreads` = as.integer(nthreads),
    display_progress = display_progress
  )
}

#' Sparse row means and variances in one pass
#'
#' Rust reimplementation of Seurat 5.6 `SparseRowMeanVar`.
#'
#' @param x Numeric non-zero values, or a `dgCMatrix`.
#' @param i Integer row indices (ignored when `x` is a `dgCMatrix`).
#' @param p Integer column pointers.
#' @param rows,cols Integer matrix dimensions.
#' @param nthreads Integer thread hint.
#' @param display_progress Logical; show progress when supported.
#'
#' @return A list with `mean` and `variance` numeric vectors.
#' @export
SparseRowMeanVar <- function(x, i, p, rows, cols, nthreads = 1L, display_progress = FALSE) {
  if (inherits(x = x, what = "sparseMatrix") && missing(i)) {
    s <- CscSlots(mat = x)
    x <- s$x
    i <- s$i
    p <- s$p
    rows <- s$nrows
    cols <- s$ncols
  }
  sparse_row_mean_var(
    x = x, i = i, p = p, nrows = rows, ncols = cols,
    `_nthreads` = as.integer(nthreads),
    `_display_progress` = display_progress
  )
}

#' Legacy sparse row means and variances
#'
#' Rust reimplementation of Seurat 5.6 `SparseRowMeanVarLegacy` (pre-5.6 VST arithmetic).
#'
#' @inheritParams SparseRowMeanVar
#'
#' @return A list with `mean` and `variance` numeric vectors.
#' @export
SparseRowMeanVarLegacy <- function(x, i, p, rows, cols, display_progress = FALSE) {
  if (inherits(x = x, what = "sparseMatrix") && missing(i)) {
    s <- CscSlots(mat = x)
    x <- s$x
    i <- s$i
    p <- s$p
    rows <- s$nrows
    cols <- s$ncols
  }
  sparse_row_mean_var_legacy(
    x = x, i = i, p = p, nrows = rows, ncols = cols,
    `_display_progress` = display_progress
  )
}

#' Clipped standardized variances for a row subset
#'
#' Rust reimplementation of Seurat 5.6 `SparseRowVarStdLegacyRows`.
#' `rows_use` is 0-based.
#'
#' @inheritParams SparseRowMeanVar
#' @param mu Numeric vector of row means.
#' @param sd Numeric vector of row standard deviations.
#' @param vmax Maximum standardized value.
#' @param rows_use Integer vector of 0-based row indices.
#'
#' @return Numeric vector of clipped row variances for `rows_use`.
#' @export
SparseRowVarStdLegacyRows <- function(
    x, i, p, mu, sd, vmax, rows, cols, rows_use,
    nthreads = 1L, display_progress = FALSE) {
  if (inherits(x = x, what = "sparseMatrix") && missing(i)) {
    s <- CscSlots(mat = x)
    x <- s$x
    i <- s$i
    p <- s$p
    rows <- s$nrows
    cols <- s$ncols
  }
  sparse_row_var_std_legacy_rows(
    x = x, i = i, p = p,
    mu = as.double(x = mu), sd = as.double(x = sd), vmax = vmax,
    `_rows` = rows, `_cols` = cols,
    rows_use = as.integer(x = rows_use),
    `_nthreads` = as.integer(nthreads),
    `_display_progress` = display_progress
  )
}

#' Bulk sparse fold-change statistics for FindAllMarkers
#'
#' Rust reimplementation of Seurat 5.6 `FindAllMarkersSparseFoldChangeStats`.
#' `groups` is 1-based.
#'
#' @inheritParams SparseRowMeanVar
#' @param groups Integer vector of 1-based group ids, length `cols`.
#' @param n_groups Number of groups.
#' @param log_normalize Logical; apply `expm1` to values before summing.
#'
#' @return A list with `group_sum`, `rest_sum`, `detected`, and `total_detected`.
#' @export
FindAllMarkersSparseFoldChangeStats <- function(
    x, i, p, rows, cols, groups, n_groups, log_normalize, nthreads = 1L) {
  if (inherits(x = x, what = "sparseMatrix") && missing(i)) {
    s <- CscSlots(mat = x)
    groups <- i
    n_groups <- p
    log_normalize <- rows
    nthreads <- if (!missing(cols)) cols else 1L
    x <- s$x
    i <- s$i
    p <- s$p
    rows <- s$nrows
    cols <- s$ncols
  }
  find_all_markers_sparse_fold_change_stats(
    x = x, i = i, p = p, rows = rows, cols = cols,
    groups = as.integer(x = groups),
    n_groups = as.integer(n_groups),
    log_normalize = log_normalize,
    `_nthreads` = as.integer(nthreads)
  )
}

#' Implicit CCA matvec without a cell-by-cell cross-product
#'
#' Rust reimplementation of Seurat 5.6 `CcaCrossprodMultiply`.
#' Computes `t(left) %*% (right %*% x)`.
#'
#' @param left,right Numeric matrices with the same number of rows.
#' @param x Numeric vector whose length matches `ncol(right)`.
#'
#' @return Numeric vector of length `ncol(left)`.
#' @export
CcaCrossprodMultiply <- function(left, right, x) {
  if (!(is.matrix(left) && is.double(left))) {
    left <- AsDoubleMatrix(x = left)
  }
  if (!(is.matrix(right) && is.double(right))) {
    right <- AsDoubleMatrix(x = right)
  }
  if (!(is.double(x) && is.null(dim(x)))) {
    x <- as.double(x = x)
  }
  cca_crossprod_multiply(left = left, right = right, x = x)
}

#' Count unique shared neighbors of integration anchors
#'
#' Rust reimplementation of Seurat 5.6 `CountAnchorSharedNeighbors`.
#' Neighbor indices are 1-based; dataset B ids are shifted by `offset`.
#'
#' @param indices_aa,indices_ab,indices_ba,indices_bb Integer/numeric neighbor index matrices.
#' @param anchor_cell1,anchor_cell2 Integer 1-based anchor cell ids.
#' @param offset Integer offset applied to dataset B indices.
#' @param k_score Integer number of neighbors to score.
#' @param nthreads Integer thread hint.
#'
#' @return Integer vector of shared-neighbor counts, one per anchor.
#' @export
CountAnchorSharedNeighbors <- function(
    indices_aa, indices_ab, indices_ba, indices_bb,
    anchor_cell1, anchor_cell2, offset, k_score, nthreads = 1L) {
  count_anchor_shared_neighbors(
    indices_aa = AsDoubleMatrix(x = indices_aa),
    indices_ab = AsDoubleMatrix(x = indices_ab),
    indices_ba = AsDoubleMatrix(x = indices_ba),
    indices_bb = AsDoubleMatrix(x = indices_bb),
    anchor_cell1 = as.integer(x = anchor_cell1),
    anchor_cell2 = as.integer(x = anchor_cell2),
    offset = as.integer(offset),
    k_score = as.integer(k_score),
    `_nthreads` = as.integer(nthreads)
  )
}

#' Approximate PCA from the Gram matrix X X'
#'
#' Rust reimplementation of Seurat 5.6 `EigenGramPCA`.
#'
#' @param object Numeric matrix of scaled data (features x cells).
#' @param npcs Integer number of principal components.
#' @param weight_by_var Logical; if `FALSE`, embeddings are divided by singular values.
#' @param nthreads Integer thread hint.
#'
#' @return A list with `loadings`, `embeddings`, and `sdev`.
#' @export
EigenGramPCA <- function(object, npcs, weight_by_var = TRUE, nthreads = 1L) {
  eigen_gram_pca(
    object = AsDoubleMatrix(x = object),
    npcs = as.integer(npcs),
    weight_by_var = weight_by_var,
    nthreads = as.integer(nthreads)
  )
}

#' SCTransform residual statistics and optional corrected counts
#'
#' Rust reimplementation of Seurat 5.6 `SCTResidualStatsAndCorrected`.
#'
#' @inheritParams SparseRowMeanVar
#' @param theta,intercept,slope Numeric per-feature model parameters.
#' @param log_umi Numeric per-cell log UMI values.
#' @param target_log_umi Numeric target log UMI for corrected counts.
#' @param min_var Minimum model variance.
#' @param residual_clip_min,residual_clip_max Residual clip range.
#' @param n_threads Integer thread hint.
#' @param compute_corrected Logical; also return corrected UMI counts.
#'
#' @return A list with `residual_mean`, `residual_variance`, and optionally `corrected`.
#' @export
SCTResidualStatsAndCorrected <- function(
    x, i, p, rows, cols, theta, intercept, slope, log_umi,
    target_log_umi, min_var, residual_clip_min, residual_clip_max,
    n_threads = 1L, compute_corrected = TRUE) {
  if (inherits(x = x, what = "sparseMatrix") && missing(i)) {
    s <- CscSlots(mat = x)
    x <- s$x
    i <- s$i
    p <- s$p
    rows <- s$nrows
    cols <- s$ncols
  }
  sct_residual_stats_and_corrected(
    x = x, i = i, p = p, rows = rows, cols = cols,
    theta = as.double(x = theta),
    intercept = as.double(x = intercept),
    slope = as.double(x = slope),
    log_umi = as.double(x = log_umi),
    target_log_umi = target_log_umi,
    min_var = min_var,
    residual_clip_min = residual_clip_min,
    residual_clip_max = residual_clip_max,
    `_n_threads` = as.integer(n_threads),
    compute_corrected = compute_corrected
  )
}

#' SCTransform Pearson residual matrix for selected features
#'
#' Rust reimplementation of Seurat 5.6 `SCTPearsonResidualMatrix`.
#' `feature_index` is 0-based.
#'
#' @inheritParams SCTResidualStatsAndCorrected
#' @param feature_index Integer 0-based feature indices to emit.
#' @param min_var Numeric scalar or per-selected-feature minimum variance.
#' @param clip_min,clip_max Residual clip range.
#' @param do_center Logical; subtract per-feature means after filling.
#' @param display_progress Logical; show progress when supported.
#'
#' @return Dense numeric matrix (selected features x cells).
#' @export
SCTPearsonResidualMatrix <- function(
    x, i, p, rows, cols, theta, intercept, slope, log_umi,
    feature_index, min_var, clip_min, clip_max, do_center = TRUE,
    n_threads = 1L, display_progress = FALSE) {
  if (inherits(x = x, what = "sparseMatrix") && missing(i)) {
    s <- CscSlots(mat = x)
    x <- s$x
    i <- s$i
    p <- s$p
    rows <- s$nrows
    cols <- s$ncols
  }
  sct_pearson_residual_matrix(
    x = x, i = i, p = p, rows = rows, cols = cols,
    theta = as.double(x = theta),
    intercept = as.double(x = intercept),
    slope = as.double(x = slope),
    log_umi = as.double(x = log_umi),
    feature_index = as.integer(x = feature_index),
    min_var = as.double(x = min_var),
    clip_min = clip_min,
    clip_max = clip_max,
    do_center = do_center,
    `_n_threads` = as.integer(n_threads),
    `_display_progress` = display_progress
  )
}

# Documented re-exports of low-level extendr entry points used by benchmarks.

#' Distances between cells and their listed neighbors
#'
#' For each row `i` of `x`, computes the Euclidean distance to the rows of `y`
#' listed in `n[[i]]`.
#'
#' @param x,y Numeric matrices of embeddings with the same number of columns
#'   (cells in rows).
#' @param n A list with one element per row of `x`, each a vector of 1-based
#'   row indices into `y`. Integer and double vectors are both accepted.
#'
#' @return A list the same length (and with the same names) as `n`, holding
#'   the distances in the same order as the neighbor indices. An empty list
#'   if `length(n) != nrow(x)`, matching Seurat.
#' @export
fast_dist <- function(x, y, n) {
  .Call(wrap__fast_dist, AsDoubleMatrix(x = x), AsDoubleMatrix(x = y), as.list(x = n))
}

#' Row sums for dgCMatrix compressed storage
#'
#' @param x Numeric non-zero values.
#' @param i Integer row indices.
#' @param rows,cols Integer matrix dimensions.
#'
#' @return Numeric vector of row sums.
#' @export
row_sum_dgcmatrix <- function(x, i, rows, cols) {
  .Call(wrap__row_sum_dgcmatrix, x, i, rows, cols)
}

#' Row means for dgCMatrix compressed storage
#'
#' @param x Numeric non-zero values.
#' @param i Integer row indices.
#' @param rows,cols Integer matrix dimensions.
#'
#' @return Numeric vector of row means.
#' @export
row_mean_dgcmatrix <- function(x, i, rows, cols) {
  .Call(wrap__row_mean_dgcmatrix, x, i, rows, cols)
}

#' Row variances for dgCMatrix compressed storage
#'
#' @param x Numeric non-zero values.
#' @param i Integer row indices.
#' @param rows,cols Integer matrix dimensions.
#'
#' @return Numeric vector of row variances.
#' @export
row_var_dgcmatrix <- function(x, i, rows, cols) {
  .Call(wrap__row_var_dgcmatrix, x, i, rows, cols)
}
