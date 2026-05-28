# High-level R wrappers matching Seurat's RcppExports API (Rust backend).

#' Extract dgCMatrix slots for extendr calls
#' @keywords internal
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
CscFromList <- function(slots) {
  methods::new(
    Class = "dgCMatrix",
    x = slots$x,
    i = slots$i,
    p = slots$p,
    Dim = as.integer(slots$Dim)
  )
}

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

#' @export
Standardize <- function(mat, display_progress = TRUE) {
  standardize(mat = mat, display_progress = display_progress)
}

#' @export
FastCov <- function(mat, center = TRUE) {
  fast_cov(mat = mat, center = center)
}

#' @export
FastCovMats <- function(mat1, mat2, center = TRUE) {
  fast_cov_mats(mat1 = mat1, mat2 = mat2, center = center)
}

#' @export
FastRBind <- function(mat1, mat2) {
  fast_rbind(mat1 = mat1, mat2 = mat2)
}

#' @export
RowVar <- function(x) {
  row_var(mat = x)
}

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

#' @export
GraphToNeighborHelper <- function(mat) {
  s <- CscSlots(mat = mat)
  graph_to_neighbor_helper(
    x = s$x, i = s$i, p = s$p, nrows = s$nrows, ncols = s$ncols
  )
}

#' @export
FastExpMean <- function(mat, display_progress) {
  s <- CscSlots(mat = mat)
  fast_exp_mean(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

#' @export
SparseRowVar <- function(mat, display_progress) {
  s <- CscSlots(mat = mat)
  sparse_row_var(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

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

#' @export
FastLogVMR <- function(mat, display_progress) {
  s <- CscSlots(mat = mat)
  fast_log_vmr(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    display_progress = display_progress
  )
}

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

#' @export
ComputeSNN <- function(nn_ranked, prune) {
  compute_snn(nn_ranked = nn_ranked, prune = prune)
}

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

#' @export
DirectSNNToFile <- function(nn_ranked, prune, display_progress, filename) {
  direct_snn_to_file(
    nn_ranked = nn_ranked,
    prune = prune,
    display_progress = display_progress,
    filename = filename
  )
}

#' @export
SNN_SmallestNonzero_Dist <- function(snn, mat, n, nearest_dist) {
  s <- CscSlots(snn)
  snn_smallest_nonzero_dist(
    x = s$x, i = s$i, p = s$p,
    nrows = s$nrows, ncols = s$ncols,
    mat = mat, n = n, nearest_dist = nearest_dist
  )
}

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
