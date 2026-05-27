# Helpers for calling Rust sparse-matrix ports from R (parity testing / migration).
# C++ remains the default via RcppExports.R.

#' Extract dgCMatrix slots for Rust extendr calls
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

#' Extract dgRMatrix slots for Rust extendr calls
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

#' Reconstruct a dgCMatrix from a Rust extendr slot list
#' @keywords internal
#' @noRd
CscFromList <- function(slots) {
  Matrix::sparseMatrix(
    i = slots$i,
    p = slots$p,
    x = slots$x,
    dims = slots$Dim,
    index1 = FALSE
  )
}

#' Rust LogNorm wrapper (parallel to `LogNorm()`)
#' @keywords internal
#' @noRd
LogNormRust <- function(data, scale_factor, display_progress = TRUE) {
  s <- CscSlots(mat = data)
  CscFromList(log_norm_rust(
    x = s$x,
    i = s$i,
    p = s$p,
    nrows = s$nrows,
    ncols = s$ncols,
    scale_factor = scale_factor,
    display_progress = display_progress
  ))
}

#' Rust RowMergeMatrices wrapper (parallel to `RowMergeMatrices()`)
#' @keywords internal
#' @noRd
RowMergeMatricesRust <- function(mat1, mat2, mat1_rownames, mat2_rownames, all_rownames) {
  s1 <- CsrSlots(mat = mat1)
  s2 <- CsrSlots(mat = mat2)
  CscFromList(row_merge_matrices_rust(
    x1 = s1$x, j1 = s1$j, p1 = s1$p, nrows1 = s1$nrows, ncols1 = s1$ncols,
    x2 = s2$x, j2 = s2$j, p2 = s2$p, nrows2 = s2$nrows, ncols2 = s2$ncols,
    mat1_rownames = mat1_rownames,
    mat2_rownames = mat2_rownames,
    all_rownames = all_rownames
  ))
}

#' Rust ReplaceColsC wrapper (parallel to `ReplaceColsC()`)
#' @keywords internal
#' @noRd
ReplaceColsRust <- function(mat, col_idx, replacement) {
  s <- CscSlots(mat = mat)
  r <- CscSlots(mat = replacement)
  CscFromList(replace_cols_rust(
    x = s$x, i = s$i, p = s$p, nrows = s$nrows, ncols = s$ncols,
    col_idx = col_idx,
    rx = r$x, ri = r$i, rp = r$p, rnrows = r$nrows, rncols = r$ncols
  ))
}

#' Rust GraphToNeighborHelper wrapper
#' @keywords internal
#' @noRd
GraphToNeighborHelperRust <- function(mat) {
  s <- CscSlots(mat = mat)
  graph_to_neighbor_helper_rust(
    x = s$x, i = s$i, p = s$p, nrows = s$nrows, ncols = s$ncols
  )
}

#' Rust FastExpMean wrapper
#' @keywords internal
#' @noRd
FastExpMeanRust <- function(mat, display_progress = TRUE) {
  s <- CscSlots(mat = mat)
  fast_exp_mean_rust(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

#' Rust SparseRowVar wrapper
#' @keywords internal
#' @noRd
SparseRowVarRust <- function(mat, display_progress = TRUE) {
  s <- CscSlots(mat = mat)
  sparse_row_var_rust(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

#' Rust FastSparseRowScale wrapper
#' @keywords internal
#' @noRd
FastSparseRowScaleRust <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  s <- CscSlots(mat = mat)
  fast_sparse_row_scale_rust(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    scale = scale, center = center, scale_max = scale_max,
    display_progress = display_progress
  )
}

#' Rust ComputeSNN wrapper
#' @keywords internal
#' @noRd
ComputeSNNRust <- function(nn_ranked, prune) {
  CscFromList(compute_snn_rust(nn_ranked = nn_ranked, prune = prune))
}

#' Rust IntegrateDataC wrapper
#' @keywords internal
#' @noRd
IntegrateDataRust <- function(integration_matrix, weights, expression_cells2) {
  im <- CscSlots(integration_matrix)
  w <- CscSlots(weights)
  ex <- CscSlots(expression_cells2)
  CscFromList(integrate_data_rust(
    ix = im$x, ii = im$i, ip = im$p, inrows = im$nrows, incols = im$ncols,
    wx = w$x, wi = w$i, wp = w$p, wnrows = w$nrows, wncols = w$ncols,
    ex = ex$x, ei = ex$i, ep = ex$p, enrows = ex$nrows, encols = ex$ncols
  ))
}

#' Rust FindWeightsC wrapper
#' @keywords internal
#' @noRd
FindWeightsRust <- function(cells2, distances, anchor_cells2, integration_matrix_rownames,
                            cell_index, anchor_score, min_dist, sd, display_progress) {
  CscFromList(find_weights_rust(
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

#' Rust fast_dist wrapper
#' @keywords internal
#' @noRd
FastDistRust <- function(x, y, n) {
  fast_dist_rust(x = x, y = y, n = n)
}

#' Rust SNN_SmallestNonzero_Dist wrapper
#' @keywords internal
#' @noRd
SNNWidthRust <- function(snn, mat, n, nearest_dist) {
  s <- CscSlots(snn)
  snn_smallest_nonzero_dist_rust(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    mat = mat, n = n, nearest_dist = nearest_dist
  )
}

#' Rust ScoreHelper wrapper
#' @keywords internal
#' @noRd
ScoreHelperRust <- function(snn, query_pca, query_dists, corrected_nns, k_snn,
                            subtract_first_nn, display_progress) {
  s <- CscSlots(snn)
  score_helper_rust(
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

#' Rust RunModularityClusteringCpp wrapper
#' @keywords internal
#' @noRd
RunModularityClusteringRust <- function(
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
  run_modularity_clustering_rust(
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
