# Shared utilities for the scRNA-seq C++ vs Rust workflow examples.
# Sourced by examples/scrna_workflow_{cpp,rust}.R and compare_scrna_workflows.R.

#' Locate repo root from examples/ or repo root cwd.
find_repo_root <- function() {
  candidates <- c(
    normalizePath(".", winslash = "/"),
    normalizePath("..", winslash = "/"),
    normalizePath("../..", winslash = "/")
  )
  for (path in candidates) {
    if (file.exists(file.path(path, "DESCRIPTION"))) {
      return(path)
    }
  }
  Sys.getenv("SEURAT_PKG_ROOT", unset = normalizePath(".", winslash = "/"))
}

bootstrap_example_env <- function() {
  root <- find_repo_root()
  setwd(root)
  bootstrap <- file.path(root, "docker/scripts/bootstrap-dev-env.R")
  if (file.exists(bootstrap)) {
    source(bootstrap, local = TRUE)
  }
  suppressPackageStartupMessages({
    devtools::load_all(recompile = FALSE, quiet = TRUE)
    library(SeuratObject)
    library(Seurat)
    library(Matrix)
    library(ggplot2)
  })
  if (!requireNamespace("seurust", quietly = TRUE)) {
    stop("seurust is not installed. Run docker/scripts/bootstrap-dev-env.R first.")
  }
  invisible(root)
}

make_backend <- function(engine = c("cpp", "rust")) {
  engine <- match.arg(engine)
  if (identical(engine, "cpp")) {
    list(
      name = "Seurat (C++/Rcpp)",
      LogNorm = function(data, scale_factor, display_progress) {
        Seurat:::LogNorm(data, scale_factor, display_progress)
      },
      row_sum_dgcmatrix = function(x, i, rows, cols) {
        Seurat:::row_sum_dgcmatrix(x, i, rows, cols)
      },
      FastExpMean = function(mat, display_progress) {
        Seurat:::FastExpMean(mat, display_progress)
      },
      FastLogVMR = function(mat, display_progress) {
        Seurat:::FastLogVMR(mat, display_progress)
      },
      SparseRowVar2 = function(mat, mu, display_progress) {
        Seurat:::SparseRowVar2(mat, mu, display_progress)
      },
      SparseRowVarStd = function(mat, mu, sd, vmax, display_progress) {
        Seurat:::SparseRowVarStd(mat, mu, sd, vmax, display_progress)
      },
      FastSparseRowScale = function(mat, scale, center, scale_max, display_progress) {
        Seurat:::FastSparseRowScale(mat, scale, center, scale_max, display_progress)
      },
      ComputeSNN = function(nn_ranked, prune) {
        Seurat:::ComputeSNN(nn_ranked, prune)
      },
      RunModularityClusteringCpp = function(SNN, modularityFunction, resolution,
                                            algorithm, nRandomStarts, nIterations,
                                            randomSeed, printOutput, edgefilename) {
        Seurat:::RunModularityClusteringCpp(
          SNN, modularityFunction, resolution, algorithm,
          nRandomStarts, nIterations, randomSeed, printOutput, edgefilename
        )
      },
      FindWeightsC = function(cells2, distances, anchor_cells2,
                              integration_matrix_rownames, cell_index,
                              anchor_score, min_dist, sd, display_progress) {
        Seurat:::FindWeightsC(
          cells2, distances, anchor_cells2, integration_matrix_rownames,
          cell_index, anchor_score, min_dist, sd, display_progress
        )
      },
      IntegrateDataC = function(integration_matrix, weights, expression_cells2) {
        Seurat:::IntegrateDataC(integration_matrix, weights, expression_cells2)
      }
    )
  } else {
    list(
      name = "seurust (Rust/extendr)",
      LogNorm = function(data, scale_factor, display_progress) {
        seurust::LogNorm(data, scale_factor, display_progress)
      },
      row_sum_dgcmatrix = function(x, i, rows, cols) {
        seurust::row_sum_dgcmatrix(x, i, rows, cols)
      },
      FastExpMean = function(mat, display_progress) {
        seurust::FastExpMean(mat, display_progress)
      },
      FastLogVMR = function(mat, display_progress) {
        seurust::FastLogVMR(mat, display_progress)
      },
      SparseRowVar2 = function(mat, mu, display_progress) {
        seurust::SparseRowVar2(mat, mu, display_progress)
      },
      SparseRowVarStd = function(mat, mu, sd, vmax, display_progress) {
        seurust::SparseRowVarStd(mat, mu, sd, vmax, display_progress)
      },
      FastSparseRowScale = function(mat, scale, center, scale_max, display_progress) {
        seurust::FastSparseRowScale(mat, scale, center, scale_max, display_progress)
      },
      ComputeSNN = function(nn_ranked, prune) {
        seurust::ComputeSNN(nn_ranked, prune)
      },
      RunModularityClusteringCpp = function(SNN, modularityFunction, resolution,
                                            algorithm, nRandomStarts, nIterations,
                                            randomSeed, printOutput, edgefilename) {
        seurust::RunModularityClusteringCpp(
          SNN, modularityFunction, resolution, algorithm,
          nRandomStarts, nIterations, randomSeed, printOutput, edgefilename
        )
      },
      FindWeightsC = function(cells2, distances, anchor_cells2,
                              integration_matrix_rownames, cell_index,
                              anchor_score, min_dist, sd, display_progress) {
        seurust::FindWeightsC(
          cells2, distances, anchor_cells2, integration_matrix_rownames,
          cell_index, anchor_score, min_dist, sd, display_progress
        )
      },
      IntegrateDataC = function(integration_matrix, weights, expression_cells2) {
        seurust::IntegrateDataC(integration_matrix, weights, expression_cells2)
      }
    )
  }
}

#' Simulate a sparse UMI count matrix resembling PBMC 10x data.
simulate_scrna_counts <- function(n_cells = 2500L, n_genes = 2000L, seed = 42L) {
  set.seed(seed)
  gene_means <- exp(rnorm(n_genes, mean = log(0.8), sd = 1.0))
  nnz <- as.integer(n_cells * n_genes * 0.12)
  i <- sample.int(n_genes, nnz, replace = TRUE)
  j <- sample.int(n_cells, nnz, replace = TRUE)
  x <- stats::rpois(nnz, lambda = pmax(gene_means[i] * 4, 0.05))
  keep <- x > 0L
  i <- i[keep]
  j <- j[keep]
  x <- x[keep]
  # Collapse duplicate (gene, cell) entries before building the dgCMatrix.
  trip <- data.frame(i = i, j = j, x = as.numeric(x))
  trip <- stats::aggregate(x ~ i + j, data = trip, FUN = sum)
  counts <- sparseMatrix(
    i = trip$i,
    j = trip$j,
    x = trip$x,
    dims = c(n_genes, n_cells),
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )
  as(counts, "dgCMatrix")
}

timed_step <- function(label, expr, timings) {
  gc(verbose = FALSE)
  t <- system.time(result <- force(expr))
  timings[[label]] <- unname(t["elapsed"])
  list(result = result, timings = timings)
}

digest_matrix <- function(mat, n = 5000L) {
  if (inherits(mat, "Matrix")) {
    v <- as.matrix(mat[seq_len(min(nrow(mat), 50L)), seq_len(min(ncol(mat), 50L))])
  } else {
    v <- mat[seq_len(min(nrow(mat), 50L)), seq_len(min(ncol(mat), 50L))]
  }
  digest::digest(v, algo = "xxhash64")
}

digest_vector <- function(x) {
  digest::digest(x, algo = "xxhash64")
}

# Native kernel steps (Rust/C++ backends only; excludes load, PCA, UMAP).
native_kernel_steps <- c(
  "02_qc_native_stats",
  "03_log_normalize",
  "04_variable_features",
  "05_scale_hvgs",
  "07_snn_graph",
  "08_clustering",
  "10_batch_integration"
)

print_timing_table <- function(timings, backend_name) {
  total <- sum(unlist(timings))
  native_total <- sum(unlist(timings[native_kernel_steps]))
  cat("\nTiming summary (", backend_name, ")\n", sep = "")
  cat(sprintf("%-28s %8s\n", "Step", "Seconds"))
  cat(strrep("-", 40), "\n", sep = "")
  for (nm in names(timings)) {
    cat(sprintf("%-28s %8.3f\n", nm, timings[[nm]]))
  }
  cat(strrep("-", 40), "\n", sep = "")
  cat(sprintf("%-28s %8.3f\n", "Total (all steps)", total))
  cat(sprintf("%-28s %8.3f\n", "Total (native kernels)", native_total))
  invisible(list(all = total, native = native_total))
}

run_scrna_workflow <- function(backend, output_file = NULL) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    install.packages("digest", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  timings <- list()
  k_param <- 20L
  dims_use <- 1:30

  # --- Load / simulate data -------------------------------------------------
  step <- timed_step("01_load_data", {
    counts <- simulate_scrna_counts(n_cells = 2500L, n_genes = 2000L, seed = 42L)
    batch <- rep(c("BatchA", "BatchB"), length.out = ncol(counts))
    seu <- CreateSeuratObject(
      counts = counts,
      project = "PBMC_demo",
      min.cells = 3,
      min.features = 50
    )
    if (ncol(seu) < 100L || nrow(seu) < 500L) {
      stop("Simulated dataset too small after CreateSeuratObject filtering.")
    }
    seu$batch <- batch[match(Cells(seu), colnames(counts))]
    seu
  }, timings)
  seu <- step$result
  timings <- step$timings

  # --- QC with native row sums ------------------------------------------------
  step <- timed_step("02_qc_native_stats", {
    mat <- GetAssayData(seu, layer = "counts")
    x <- slot(mat, "x")
    i <- slot(mat, "i")
    n_count <- backend$row_sum_dgcmatrix(x, i, nrow(mat), ncol(mat))
    seu$nCount_RNA <- n_count[colnames(seu)]
    seu$nFeature_RNA <- Matrix::colSums(mat > 0)
    seu <- subset(seu, subset = nFeature_RNA > 150 & nFeature_RNA < 2200)
    seu
  }, timings)
  seu <- step$result
  timings <- step$timings

  # --- Log normalization (native) ---------------------------------------------
  step <- timed_step("03_log_normalize", {
    mat <- GetAssayData(seu, layer = "counts")
    norm <- backend$LogNorm(mat, scale_factor = 1e4, display_progress = FALSE)
    rownames(norm) <- rownames(mat)
    colnames(norm) <- colnames(mat)
    seu <- SetAssayData(seu, layer = "data", new.data = norm)
    seu
  }, timings)
  seu <- step$result
  timings <- step$timings

  # --- Variable features (native VST-style) -----------------------------------
  step <- timed_step("04_variable_features", {
    mat <- GetAssayData(seu, layer = "data")
    hvf <- data.frame(mean = Matrix::rowMeans(mat))
    hvf$variance <- backend$SparseRowVar2(mat, mu = hvf$mean, display_progress = FALSE)
    clip.max <- sqrt(ncol(mat))
    fit <- stats::loess(
      log10(variance) ~ log10(mean),
      data = hvf[hvf$variance > 0, , drop = FALSE],
      span = 0.3
    )
    hvf$variance.expected <- 0
    hvf$variance.expected[hvf$variance > 0] <- 10^fit$fitted
    hvf$variance.standardized <- backend$SparseRowVarStd(
      mat,
      mu = hvf$mean,
      sd = sqrt(hvf$variance.expected),
      vmax = clip.max,
      display_progress = FALSE
    )
    top_hvf <- rownames(hvf)[order(hvf$variance.standardized, decreasing = TRUE)][1:2000]
    VariableFeatures(seu) <- top_hvf
    seu
  }, timings)
  seu <- step$result
  timings <- step$timings
  hvf <- VariableFeatures(seu)

  # --- Scale HVGs (native) ----------------------------------------------------
  step <- timed_step("05_scale_hvgs", {
    mat <- GetAssayData(seu, layer = "data")[hvf, , drop = FALSE]
    scaled <- backend$FastSparseRowScale(
      mat,
      scale = TRUE,
      center = TRUE,
      scale_max = 10,
      display_progress = FALSE
    )
    scaled_dense <- as.matrix(scaled)
    rownames(scaled_dense) <- rownames(mat)
    colnames(scaled_dense) <- colnames(mat)
    seu <- SetAssayData(seu, layer = "scale.data", new.data = scaled_dense)
    seu
  }, timings)
  seu <- step$result
  timings <- step$timings

  # --- PCA (Seurat + irlba; same backend for both scripts) --------------------
  step <- timed_step("06_pca", {
    RunPCA(
      seu,
      features = hvf,
      npcs = 30,
      verbose = FALSE,
      rev.pca = TRUE
    )
  }, timings)
  seu <- step$result
  timings <- step$timings

  # --- Neighbors: shared kNN index, native SNN graph --------------------------
  step <- timed_step("07_snn_graph", {
    emb <- Embeddings(seu, reduction = "pca")[, dims_use, drop = FALSE]
    nn <- Seurat:::NNHelper(
      data = emb,
      k = k_param,
      method = "annoy",
      n.trees = 50,
      searchtype = "standard",
      metric = "euclidean"
    )
    nn_ranked <- Indices(nn)
    snn <- backend$ComputeSNN(nn_ranked = nn_ranked, prune = 1 / 15)
    rownames(snn) <- colnames(snn) <- rownames(emb)
    list(snn = snn, nn_ranked = nn_ranked)
  }, timings)
  snn <- step$result$snn
  nn_ranked <- step$result$nn_ranked
  timings <- step$timings

  # --- Louvain clustering (native modularity optimizer) -------------------------
  step <- timed_step("08_clustering", {
    backend$RunModularityClusteringCpp(
      SNN = snn,
      modularityFunction = 1L,
      resolution = 0.8,
      algorithm = 1L,
      nRandomStarts = 10L,
      nIterations = 10L,
      randomSeed = 42L,
      printOutput = FALSE,
      edgefilename = ""
    )
  }, timings)
  clusters <- step$result
  timings <- step$timings
  names(clusters) <- rownames(snn)
  seu$seurat_clusters <- factor(clusters)

  # --- UMAP for visualization (shared Seurat path) ----------------------------
  step <- timed_step("09_umap", {
    seu <- RunUMAP(seu, dims = dims_use, verbose = FALSE, seed.use = 42L)
    seu
  }, timings)
  seu <- step$result
  timings <- step$timings

  # --- Mini batch-integration demo (native FindWeightsC + IntegrateDataC) -----
  step <- timed_step("10_batch_integration", {
    # Split by batch (typical multi-sample workflow); native integration kernels
    # use the same small fixtures as tests/testthat/test_rust_cpp_parity_snn_integration.R.
    obj.list <- SplitObject(seu, split.by = "batch")
    invisible(list(ref = obj.list[[1]], qry = obj.list[[2]]))

    cell_index <- matrix(c(1, 2, 2, 1), nrow = 2, byrow = TRUE)
    storage.mode(cell_index) <- "double"
    weights <- backend$FindWeightsC(
      cells2 = as.double(c(0, 1)),
      distances = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, byrow = TRUE),
      anchor_cells2 = c("Gene1", "Gene2"),
      integration_matrix_rownames = c("Gene1", "Gene2", "Gene1"),
      cell_index = cell_index,
      anchor_score = c(1, 0.5, 0.8),
      min_dist = 0,
      sd = 1,
      display_progress = FALSE
    )
    expr <- as(
      sparseMatrix(
        i = c(1, 2, 3, 1, 2),
        j = c(1, 1, 1, 2, 2),
        x = c(1, 2, 3, 4, 5),
        dims = c(3L, 2L)
      ),
      "dgCMatrix"
    )
    im <- as(
      sparseMatrix(
        i = c(1, 2, 1),
        j = c(1, 2, 1),
        x = c(0.5, 0.3, 0.2),
        dims = c(2L, 2L)
      ),
      "dgCMatrix"
    )
    w <- as(
      sparseMatrix(
        i = c(1, 2, 1),
        j = c(1, 2, 3),
        x = c(0.4, 0.6, 0.1),
        dims = c(2L, 3L)
      ),
      "dgCMatrix"
    )
    integrated <- backend$IntegrateDataC(
      integration_matrix = im,
      weights = w,
      expression_cells2 = expr
    )
    list(
      anchor_count = 2L,
      weights_digest = digest_matrix(weights),
      integrated_digest = digest_matrix(integrated)
    )
  }, timings)
  integration_summary <- step$result
  timings <- step$timings

  # --- Summary outputs --------------------------------------------------------
  umap_plot <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste("Clusters (", backend$name, ")", sep = "")) +
    theme(legend.position = "none")

  results <- list(
    backend = backend$name,
    n_cells = ncol(seu),
    n_clusters = length(unique(clusters)),
    cluster_table = as.integer(table(clusters)),
    cluster_digest = digest_vector(clusters),
    norm_digest = digest_matrix(GetAssayData(seu, layer = "data")),
    snn_digest = digest_matrix(snn),
    snn_nnz = length(slot(snn, "x")),
    umap_digest = digest_matrix(Embeddings(seu, "umap")),
    integration = integration_summary,
    timings = timings,
    total_native_seconds = sum(unlist(timings)),
    total_kernel_seconds = sum(unlist(timings[native_kernel_steps]))
  )

  print_timing_table(timings, backend$name)
  cat("\nWorkflow outputs:\n")
  cat("  Cells after QC:", results$n_cells, "\n")
  cat("  Clusters found:", results$n_clusters, "\n")
  cat("  Cluster sizes:", paste(results$cluster_table, collapse = ", "), "\n")
  cat("  SNN non-zero entries:", results$snn_nnz, "\n")
  cat("  Integration anchors:", integration_summary$anchor_count, "\n")

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(results, output_file)
    ggsave(
      sub("\\.rds$", "_umap.png", output_file),
      plot = umap_plot,
      width = 6,
      height = 5,
      dpi = 120
    )
    cat("  Saved:", output_file, "\n")
  }

  invisible(list(seurat = seu, results = results, plot = umap_plot))
}
