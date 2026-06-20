# Capture and compare upstream script outputs between Seurat C++ and RSeurat runs.

ensure_digest_pkg <- function() {
  if (!requireNamespace("digest", quietly = TRUE)) {
    install.packages("digest", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

digest_vector <- function(x) {
  ensure_digest_pkg()
  digest::digest(as.vector(x), algo = "xxhash64")
}

digest_matrix <- function(mat, n_rows = 100L, n_cols = 100L) {
  ensure_digest_pkg()
  if (is.null(mat) || length(mat) == 0L) {
    return(NA_character_)
  }
  nr <- min(nrow(mat), n_rows)
  nc <- min(ncol(mat), n_cols)
  sample <- mat[seq_len(nr), seq_len(nc), drop = FALSE]
  digest::digest(sample, algo = "xxhash64")
}

digest_rounded_matrix <- function(mat, digits = 8L, n_rows = 100L, n_cols = 100L) {
  ensure_digest_pkg()
  if (is.null(mat) || length(mat) == 0L) {
    return(NA_character_)
  }
  nr <- min(nrow(mat), n_rows)
  nc <- min(ncol(mat), n_cols)
  sample <- mat[seq_len(nr), seq_len(nc), drop = FALSE]
  digest::digest(as.vector(round(unname(sample), digits = digits)), algo = "xxhash64")
}

# irlba/SVD PC directions are defined only up to sign and scale; align before digesting.
normalize_embedding_columns <- function(mat) {
  if (is.null(mat) || length(mat) == 0L || ncol(mat) == 0L) {
    return(mat)
  }
  scaled <- scale(mat, center = FALSE, scale = sqrt(colSums(mat^2)))
  if (is.null(scaled)) {
    return(mat * 0)
  }
  unname(scaled)
}

canonicalize_embedding_matrix <- function(mat) {
  if (is.null(mat) || length(mat) == 0L || ncol(mat) == 0L) {
    return(mat)
  }
  out <- normalize_embedding_columns(mat)
  for (j in seq_len(ncol(out))) {
    col <- out[, j]
    peak <- which.max(abs(col))
    if (length(peak) == 0L || col[peak] < 0) {
      out[, j] <- -col
    }
  }
  out
}

embedding_sample <- function(mat, n_rows = 100L, n_cols = 100L) {
  if (is.null(mat) || length(mat) == 0L) {
    return(numeric())
  }
  nr <- min(nrow(mat), n_rows)
  nc <- min(ncol(mat), n_cols)
  as.vector(unname(mat[seq_len(nr), seq_len(nc), drop = FALSE]))
}

embedding_fingerprint <- function(mat, tolerance = 1e-4) {
  ensure_digest_pkg()
  values <- embedding_sample(mat)
  list(
    values = values,
    tolerance = tolerance,
    digest = digest::digest(round(values, digits = ceiling(-log10(tolerance))), algo = "xxhash64")
  )
}

digest_file <- function(path) {
  ensure_digest_pkg()
  if (!file.exists(path)) {
    return(NA_character_)
  }
  digest::digest(file = path, algo = "xxhash64")
}

digest_seurat_object <- function(obj) {
  out <- list()
  if (!inherits(obj, "Seurat")) {
    return(out)
  }

  if (length(as.character(Idents(obj))) > 0L) {
    out$cluster_digest <- digest_vector(as.character(Idents(obj)))
    out$n_clusters <- length(unique(as.character(Idents(obj))))
  }

  graph_names <- grep("_snn$", names(obj@graphs), value = TRUE)
  if (length(graph_names) > 0L) {
    snn <- obj@graphs[[graph_names[[1]]]]
    out$snn_digest <- digest_matrix(snn)
    out$snn_nnz <- length(snn@x)
  }

  scale_data <- tryCatch(
    GetAssayData(obj, layer = "scale.data"),
    error = function(e) NULL
  )
  if (!is.null(scale_data) && length(scale_data) > 0L) {
    out$scaled_rounded_digest <- digest_rounded_matrix(scale_data)
  }

  for (reduction in c("harmony", "umap", "pca")) {
    if (reduction %in% names(obj@reductions)) {
      embed <- Embeddings(obj, reduction = reduction)
      if (reduction %in% c("pca", "harmony")) {
        embed <- canonicalize_embedding_matrix(embed)
      }
      field <- paste0(reduction, "_digest")
      out[[field]] <- if (reduction %in% c("pca", "harmony")) {
        embedding_fingerprint(embed)
      } else {
        digest_rounded_matrix(embed)
      }
    }
  }

  out
}

SCRIPT_OUTPUTS <- list(
  "01_pbmc_satija_tutorial.R" = list(
    rds = c("../output/pbmc3k_final_002.rds"),
    csv = c("../output/pbmc3k_markers_002.csv")
  ),
  "02_gbm_seurat_adapted.R" = list(
    rds = c(
      "../output/gbm_unharmonized_002.rds",
      "../output/gbm_harmonized_002.rds"
    ),
    csv = c("../output/qc_summary_002.csv")
  ),
  "03_gbmap_exploration.R" = list(
    csv = c(
      "../output/celltype_crosstab_002.csv",
      "../output/donor_composition_002.csv"
    )
  )
)

capture_script_outputs <- function(script_name, code_dir) {
  spec <- SCRIPT_OUTPUTS[[script_name]]
  if (is.null(spec)) {
    return(list())
  }

  digests <- list()
  for (rds_rel in spec$rds %||% character()) {
    path <- normalizePath(file.path(code_dir, rds_rel), winslash = "/", mustWork = FALSE)
    label <- sub("\\.rds$", "", basename(path))
    if (!file.exists(path)) {
      digests[[paste0(label, "_missing")]] <- TRUE
      next
    }
    obj <- readRDS(path)
    obj_digests <- digest_seurat_object(obj)
    for (nm in names(obj_digests)) {
      digests[[paste0(label, "_", nm)]] <- obj_digests[[nm]]
    }
  }

  for (csv_rel in spec$csv %||% character()) {
    path <- normalizePath(file.path(code_dir, csv_rel), winslash = "/", mustWork = FALSE)
    label <- sub("\\.csv$", "", basename(path))
    digests[[paste0(label, "_file_digest")]] <- digest_file(path)
  }

  digests
}

`%||%` <- function(x, y) if (is.null(x)) y else x

compare_script_outputs <- function(script_name, cpp_digests, rust_digests) {
  keys <- sort(unique(c(names(cpp_digests), names(rust_digests))))
  if (length(keys) == 0L) {
    return(invisible(TRUE))
  }

  exact_suffixes <- c(
    "_cluster_digest",
    "_n_clusters",
    "_scaled_rounded_digest",
    "_pca_digest",
    "_harmony_digest",
    "_snn_digest",
    "_snn_nnz",
    "_file_digest"
  )
  is_exact <- function(key) {
    any(vapply(exact_suffixes, function(suffix) endsWith(key, suffix), logical(1)))
  }
  is_embedding <- function(key) {
    endsWith(key, "_pca_digest") || endsWith(key, "_harmony_digest")
  }
  compare_values <- function(key, cpp_val, rust_val) {
    if (is_embedding(key) && is.list(cpp_val) && is.list(rust_val)) {
      if (!identical(length(cpp_val$values), length(rust_val$values))) {
        return(FALSE)
      }
      if (length(cpp_val$values) == 0L) {
        return(TRUE)
      }
      max(abs(cpp_val$values - rust_val$values), na.rm = TRUE) <= cpp_val$tolerance
    } else {
      identical(cpp_val, rust_val)
    }
  }
  format_value <- function(val) {
    if (is.list(val) && !is.null(val$digest)) {
      return(paste0(val$digest, " (tol=", val$tolerance, ")"))
    }
    as.character(val)
  }

  cat("\n==> Output parity: ", script_name, "\n", sep = "")
  all_ok <- TRUE
  for (key in keys) {
    cpp_val <- cpp_digests[[key]]
    rust_val <- rust_digests[[key]]
    ok <- compare_values(key, cpp_val, rust_val)
    exact <- is_exact(key)
    if (exact) {
      all_ok <- all_ok && ok
    }
    status <- if (ok) {
      "OK"
    } else if (exact) {
      "MISMATCH"
    } else {
      "MISMATCH (informational)"
    }
    cat(sprintf("  %-40s %s\n", key, status))
    if (!ok) {
      cat("    Seurat C++:", format_value(cpp_val), "\n")
      cat("    RSeurat   :", format_value(rust_val), "\n")
    }
  }

  if (!all_ok) {
    stop(
      "Output mismatch between Seurat C++ and RSeurat for ",
      script_name,
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
