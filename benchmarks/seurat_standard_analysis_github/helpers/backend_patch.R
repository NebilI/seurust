# Swap Seurat C++ RcppExports entry points with RSeurat Rust implementations.
# Lets upstream analysis scripts run unchanged while benchmarking backends.

.patch_store <- new.env(parent = emptyenv())

PORTED_NATIVE_FUNS <- c(
  "LogNorm",
  "Standardize",
  "FastSparseRowScale",
  "FastSparseRowScaleWithKnownStats",
  "FastCov",
  "FastCovMats",
  "FastRBind",
  "RowVar",
  "RowMergeMatrices",
  "ReplaceColsC",
  "GraphToNeighborHelper",
  "FastExpMean",
  "SparseRowVar",
  "SparseRowVar2",
  "SparseRowVarStd",
  "FastLogVMR",
  "RunUMISampling",
  "RunUMISamplingPerCell",
  "ComputeSNN",
  "IntegrateDataC",
  "FindWeightsC",
  "ScoreHelper",
  "WriteEdgeFile",
  "DirectSNNToFile",
  "SNN_SmallestNonzero_Dist",
  "RunModularityClusteringCpp",
  "row_sum_dgcmatrix"
)

patch_seurat_backend <- function(engine = c("cpp", "rust")) {
  engine <- match.arg(engine)
  if (!requireNamespace("RSeurat", quietly = TRUE)) {
    stop("RSeurat is not installed.", call. = FALSE)
  }

  seurat_ns <- asNamespace("Seurat")
  rseurat_ns <- asNamespace("RSeurat")
  patched <- character()

  for (fn in PORTED_NATIVE_FUNS) {
    if (!exists(fn, envir = seurat_ns, inherits = FALSE)) {
      next
    }
    if (!exists(fn, envir = .patch_store, inherits = FALSE)) {
      assign(fn, get(fn, envir = seurat_ns), envir = .patch_store)
    }
    replacement <- if (identical(engine, "rust")) {
      if (!exists(fn, envir = rseurat_ns, inherits = FALSE)) {
        next
      }
      get(fn, envir = rseurat_ns)
    } else {
      get(fn, envir = .patch_store)
    }
    if (bindingIsLocked(fn, seurat_ns)) {
      unlockBinding(fn, seurat_ns)
    }
    assign(fn, replacement, envir = seurat_ns)
    lockBinding(fn, seurat_ns)
    patched <- c(patched, fn)
  }

  attr(patched, "engine") <- engine
  invisible(patched)
}

backend_label <- function(engine) {
  if (identical(engine, "rust")) "RSeurat Rust" else "Seurat C++"
}
