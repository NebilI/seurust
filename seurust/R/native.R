# High-level R wrappers matching Seurat's RcppExports API (Rust backend).

#' Extract dgCMatrix slots for extendr calls
#' @keywords internal
#' @noRd
CscSlots <- function(mat) {
  if (!inherits(x = mat, what = "dgCMatrix")) {
    mat <- as(object = mat, Class = "dgCMatrix")
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
    mat <- as(object = mat, Class = "RsparseMatrix")
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
Standardize <- function(mat, display_progress = TRUE) {
  standardize(mat = mat, display_progress = display_progress)
}

#' Fast covariance of a dense matrix
#'
#' @param mat Numeric matrix.
#' @param center Logical; center columns before computing covariance.
#'
#' @return A symmetric numeric covariance matrix.
#' @export
FastCov <- function(mat, center = TRUE) {
  fast_cov(mat = mat, center = center)
}

#' Fast cross-covariance of two dense matrices
#'
#' @param mat1,mat2 Numeric matrices with the same number of rows.
#' @param center Logical; center columns before computing covariance.
#'
#' @return A numeric cross-covariance matrix.
#' @export
FastCovMats <- function(mat1, mat2, center = TRUE) {
  fast_cov_mats(mat1 = mat1, mat2 = mat2, center = center)
}

#' Fast row-bind of two dense matrices
#'
#' @param mat1,mat2 Numeric matrices with the same number of columns.
#'
#' @return A numeric matrix stacking `mat1` and `mat2` by row.
#' @export
FastRBind <- function(mat1, mat2) {
  fast_rbind(mat1 = mat1, mat2 = mat2)
}

#' Row variances of a dense matrix
#'
#' @param x Numeric matrix.
#'
#' @return Numeric vector of row variances.
#' @export
RowVar <- function(x) {
  row_var(mat = x)
}

#' Merge two sparse matrices by shared row names
#'
#' @param mat1,mat2 Sparse matrices in row-compressed form (or coercible).
#' @param mat1_rownames,mat2_rownames Character vectors of row names.
#' @param all_rownames Character vector of the union of row names.
#'
#' @return A `dgCMatrix` containing the merged rows.
#' @export
RowMergeMatrices <- function(mat1, mat2, mat1_rownames, mat2_rownames, all_rownames) {
  s1 <- CsrSlots(mat = mat1)
  s2 <- CsrSlots(mat = mat2)
  CscFromList(row_merge_matrices(
    x1 = s1$x, j1 = s1$j, p1 = s1$p, nrows1 = s1$nrows, ncols1 = s1$ncols,
    x2 = s2$x, j2 = s2$j, p2 = s2$p, nrows2 = s2$nrows, ncols2 = s2$ncols,
    mat1_rownames = mat1_rownames,
    mat2_rownames = mat2_rownames,
    all_rownames = all_rownames
  ))
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
  CscFromList(replace_cols(
    x = s$x, i = s$i, p = s$p, nrows = s$nrows, ncols = s$ncols,
    col_idx = col_idx,
    rx = r$x, ri = r$i, rp = r$p, rnrows = r$nrows, rncols = r$ncols
  ))
}

#' Convert a sparse graph to neighbor index lists
#'
#' @param mat A sparse adjacency/`Graph` matrix (`dgCMatrix`).
#'
#' @return Neighbor index information matching Seurat's helper.
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
    mu = mu,
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
    mu = mu, sd = sd, vmax = vmax,
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
FastSparseRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  s <- CscSlots(mat = mat)
  fast_sparse_row_scale(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
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
    mu = mu, sigma = sigma,
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
    sample_val = sample_val, upsample = upsample,
    `_display_progress` = display_progress
  )
  CscFromList(out)
}

#' Build a shared nearest-neighbor (SNN) graph
#'
#' @param nn_ranked Integer matrix of neighbor ranks.
#' @param prune Numeric pruning threshold.
#'
#' @return A sparse SNN graph as a `dgCMatrix`.
#' @export
ComputeSNN <- function(nn_ranked, prune) {
  compute_snn(nn_ranked = nn_ranked, prune = prune)
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
  CscFromList(integrate_data(
    ix = im$x, ii = im$i, ip = im$p, inrows = im$nrows, incols = im$ncols,
    wx = w$x, wi = w$i, wp = w$p, wnrows = w$nrows, wncols = w$ncols,
    ex = ex$x, ei = ex$i, ep = ex$p, enrows = ex$nrows, encols = ex$ncols
  ))
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
                         cell_index, anchor_score, min_dist, sd, display_progress) {
  CscFromList(find_weights(
    cells2 = cells2,
    distances = distances,
    anchor_cells2 = anchor_cells2,
    integration_matrix_rownames = integration_matrix_rownames,
    cell_index = cell_index,
    anchor_score = anchor_score,
    min_dist = min_dist,
    sd = sd,
    display_progress = display_progress
  ))
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
                        subtract_first_nn, display_progress) {
  s <- CscSlots(snn)
  score_helper(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    query_pca = query_pca,
    query_dists = query_dists,
    corrected_nns = corrected_nns,
    k_snn = k_snn,
    subtract_first_nn = subtract_first_nn,
    display_progress = display_progress
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
#' @param nn_ranked Integer matrix of neighbor ranks.
#' @param prune Numeric pruning threshold.
#' @param display_progress Logical; show progress when supported.
#' @param filename Output path.
#'
#' @return Invisibly returns a status code or `NULL`, matching Seurat.
#' @export
DirectSNNToFile <- function(nn_ranked, prune, display_progress, filename) {
  direct_snn_to_file(
    nn_ranked = nn_ranked,
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
    mat = mat, n = n, nearest_dist = nearest_dist
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
    edgefilename = "") {
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
    edgefilename = edgefilename
  )
}

# Documented re-exports of low-level extendr entry points used by benchmarks.

#' Fast pairwise distances for neighbor search
#'
#' @param x,y Numeric matrices of embeddings.
#' @param n Integer number of neighbors.
#'
#' @return Neighbor indices/distances matching Seurat's `fast_dist`.
#' @export
fast_dist <- function(x, y, n) {
  .Call(wrap__fast_dist, x, y, n)
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
