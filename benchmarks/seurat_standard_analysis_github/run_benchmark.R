#!/usr/bin/env Rscript
# Benchmark ajtimon/seurat-standard-analysis PBMC tutorial with Seurat C++ vs RSeurat Rust kernels.

UPSTREAM_REPO_URL <- "https://github.com/ajtimon/seurat-standard-analysis"
UPSTREAM_SCRIPT <- "code/01_pbmc_satija_tutorial.R"
UPSTREAM_WORKFLOW_URL <- paste0(UPSTREAM_REPO_URL, "/blob/master/", UPSTREAM_SCRIPT)
PBMC3K_DATA_URL <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

find_repo_root <- function() {
  candidates <- c(
    normalizePath(".", winslash = "/", mustWork = FALSE),
    normalizePath("..", winslash = "/", mustWork = FALSE),
    normalizePath("../..", winslash = "/", mustWork = FALSE),
    normalizePath("../../..", winslash = "/", mustWork = FALSE)
  )
  for (path in candidates) {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        dir.exists(file.path(path, "RSeurat"))) {
      return(path)
    }
  }
  Sys.getenv("SEURAT_PKG_ROOT", unset = "/workspace")
}

repo_root <- find_repo_root()
setwd(repo_root)

system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")

if (!requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest", repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  devtools::load_all(recompile = FALSE, quiet = TRUE)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(RSeurat)
})

benchmark_root <- file.path(repo_root, "benchmarks/seurat_standard_analysis_github")
data_root <- Sys.getenv(
  "SEURAT_STANDARD_ANALYSIS_DATA_DIR",
  unset = file.path(benchmark_root, "data")
)
output_root <- Sys.getenv(
  "SEURAT_STANDARD_ANALYSIS_OUTPUT_DIR",
  unset = file.path(benchmark_root, "output")
)
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

pbmc_matrix_dir <- file.path(data_root, "filtered_gene_bc_matrices/hg19")
pbmc_archive <- file.path(data_root, "pbmc3k_filtered_gene_bc_matrices.tar.gz")

download_pbmc3k <- function() {
  matrix_file <- file.path(pbmc_matrix_dir, "matrix.mtx")
  gz_matrix_file <- paste0(matrix_file, ".gz")
  if (file.exists(matrix_file) || file.exists(gz_matrix_file)) {
    return(invisible(pbmc_matrix_dir))
  }

  message("==> Downloading PBMC 3K data from 10x Genomics...")
  utils::download.file(PBMC3K_DATA_URL, pbmc_archive, mode = "wb", quiet = FALSE)
  utils::untar(pbmc_archive, exdir = data_root)

  if (!dir.exists(pbmc_matrix_dir)) {
    stop("PBMC 3K archive did not contain expected directory: ", pbmc_matrix_dir)
  }
  invisible(pbmc_matrix_dir)
}

timed_step <- function(label, expr, timings) {
  gc(verbose = FALSE)
  elapsed <- system.time(result <- force(expr))[["elapsed"]]
  timings[[label]] <- unname(elapsed)
  list(result = result, timings = timings)
}

digest_matrix <- function(mat, n_rows = 100L, n_cols = 100L) {
  rows <- seq_len(min(nrow(mat), n_rows))
  cols <- seq_len(min(ncol(mat), n_cols))
  sample <- if (inherits(mat, "Matrix")) {
    as.matrix(mat[rows, cols, drop = FALSE])
  } else {
    mat[rows, cols, drop = FALSE]
  }
  digest::digest(sample, algo = "xxhash64")
}

digest_rounded_matrix <- function(mat, digits = 8L, n_rows = 100L, n_cols = 100L) {
  rows <- seq_len(min(nrow(mat), n_rows))
  cols <- seq_len(min(ncol(mat), n_cols))
  sample <- if (inherits(mat, "Matrix")) {
    as.matrix(mat[rows, cols, drop = FALSE])
  } else {
    mat[rows, cols, drop = FALSE]
  }
  digest::digest(round(sample, digits = digits), algo = "xxhash64")
}

digest_vector <- function(x) {
  digest::digest(as.vector(x), algo = "xxhash64")
}

make_backend <- function(engine = c("cpp", "rust")) {
  engine <- match.arg(engine)
  if (identical(engine, "cpp")) {
    list(
      name = "Seurat C++",
      LogNorm = Seurat:::LogNorm,
      SparseRowVar2 = Seurat:::SparseRowVar2,
      SparseRowVarStd = Seurat:::SparseRowVarStd,
      FastSparseRowScale = Seurat:::FastSparseRowScale,
      ComputeSNN = Seurat:::ComputeSNN,
      RunModularityClusteringCpp = Seurat:::RunModularityClusteringCpp
    )
  } else {
    list(
      name = "RSeurat Rust",
      LogNorm = RSeurat::LogNorm,
      SparseRowVar2 = RSeurat::SparseRowVar2,
      SparseRowVarStd = RSeurat::SparseRowVarStd,
      FastSparseRowScale = RSeurat::FastSparseRowScale,
      ComputeSNN = RSeurat::ComputeSNN,
      RunModularityClusteringCpp = RSeurat::RunModularityClusteringCpp
    )
  }
}

select_variable_features <- function(seu, backend, nfeatures = 2000L) {
  mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  hvf <- data.frame(mean = Matrix::rowMeans(mat))
  hvf$variance <- backend$SparseRowVar2(mat, mu = hvf$mean, display_progress = FALSE)

  fit_rows <- is.finite(hvf$mean) & hvf$mean > 0 & is.finite(hvf$variance) & hvf$variance > 0
  fit <- stats::loess(
    log10(variance) ~ log10(mean),
    data = hvf[fit_rows, , drop = FALSE],
    span = 0.3
  )
  hvf$variance.expected <- 0
  hvf$variance.expected[fit_rows] <- 10^stats::predict(fit, hvf[fit_rows, , drop = FALSE])
  hvf$variance.expected[!is.finite(hvf$variance.expected) | hvf$variance.expected <= 0] <- NA_real_

  hvf$variance.standardized <- backend$SparseRowVarStd(
    mat,
    mu = hvf$mean,
    sd = sqrt(hvf$variance.expected),
    vmax = sqrt(ncol(mat)),
    display_progress = FALSE
  )
  hvf$variance.standardized[!is.finite(hvf$variance.standardized)] <- -Inf

  top <- rownames(hvf)[order(hvf$variance.standardized, decreasing = TRUE)]
  top[seq_len(min(nfeatures, length(top)))]
}

run_workflow <- function(backend, counts, output_file) {
  timings <- list()
  dims_use <- 1:10
  k_param <- 20L

  step <- timed_step("01_create_and_qc", {
    seu <- CreateSeuratObject(
      counts = counts,
      project = "pbmc3k",
      min.cells = 3,
      min.features = 200
    )
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  }, timings)
  seu <- step$result
  timings <- step$timings

  step <- timed_step("02_log_normalize", {
    mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
    norm <- backend$LogNorm(mat, scale_factor = 1e4, display_progress = FALSE)
    rownames(norm) <- rownames(mat)
    colnames(norm) <- colnames(mat)
    SetAssayData(seu, assay = "RNA", layer = "data", new.data = norm)
  }, timings)
  seu <- step$result
  timings <- step$timings

  step <- timed_step("03_variable_features", {
    VariableFeatures(seu) <- select_variable_features(seu, backend)
    seu
  }, timings)
  seu <- step$result
  timings <- step$timings
  hvf <- VariableFeatures(seu)

  step <- timed_step("04_scale_data", {
    mat <- GetAssayData(seu, assay = "RNA", layer = "data")
    scaled <- backend$FastSparseRowScale(
      mat,
      scale = TRUE,
      center = TRUE,
      scale_max = 10,
      display_progress = FALSE
    )
    scaled <- as.matrix(scaled)
    rownames(scaled) <- rownames(mat)
    colnames(scaled) <- colnames(mat)
    SetAssayData(seu, assay = "RNA", layer = "scale.data", new.data = scaled)
  }, timings)
  seu <- step$result
  timings <- step$timings

  step <- timed_step("05_pca", {
    RunPCA(seu, features = hvf, npcs = 30, verbose = FALSE)
  }, timings)
  seu <- step$result
  timings <- step$timings

  step <- timed_step("06_snn_graph", {
    emb <- Embeddings(seu, reduction = "pca")[, dims_use, drop = FALSE]
    nn <- Seurat:::NNHelper(
      data = emb,
      k = k_param,
      method = "annoy",
      n.trees = 50,
      searchtype = "standard",
      metric = "euclidean"
    )
    snn <- backend$ComputeSNN(nn_ranked = Indices(nn), prune = 1 / 15)
    rownames(snn) <- colnames(snn) <- rownames(emb)
    snn
  }, timings)
  snn <- step$result
  timings <- step$timings

  step <- timed_step("07_clustering", {
    clusters <- backend$RunModularityClusteringCpp(
      SNN = snn,
      modularityFunction = 1L,
      resolution = 0.5,
      algorithm = 1L,
      nRandomStarts = 10L,
      nIterations = 10L,
      randomSeed = 42L,
      printOutput = FALSE,
      edgefilename = ""
    )
    names(clusters) <- rownames(snn)
    clusters
  }, timings)
  clusters <- step$result
  timings <- step$timings
  seu$seurat_clusters <- factor(clusters)

  step <- timed_step("08_umap", {
    RunUMAP(
      seu,
      dims = dims_use,
      verbose = FALSE,
      seed.use = 42L,
      umap.method = "uwot",
      metric = "cosine"
    )
  }, timings)
  seu <- step$result
  timings <- step$timings

  results <- list(
    upstream_repo_url = UPSTREAM_REPO_URL,
    upstream_script = UPSTREAM_SCRIPT,
    upstream_workflow_url = UPSTREAM_WORKFLOW_URL,
    backend = backend$name,
    n_cells = ncol(seu),
    n_genes = nrow(seu),
    n_clusters = length(unique(clusters)),
    cluster_table = as.integer(table(clusters)),
    variable_features_digest = digest_vector(hvf),
    normalized_digest = digest_matrix(GetAssayData(seu, assay = "RNA", layer = "data")),
    scaled_digest = digest_matrix(GetAssayData(seu, assay = "RNA", layer = "scale.data")),
    scaled_rounded_digest = digest_rounded_matrix(
      GetAssayData(seu, assay = "RNA", layer = "scale.data")
    ),
    snn_digest = digest_matrix(snn),
    snn_nnz = length(slot(snn, "x")),
    cluster_digest = digest_vector(clusters),
    umap_digest = digest_matrix(Embeddings(seu, "umap")),
    timings = timings,
    native_kernel_seconds = sum(unlist(timings[c(
      "02_log_normalize",
      "03_variable_features",
      "04_scale_data",
      "06_snn_graph",
      "07_clustering"
    )])),
    total_seconds = sum(unlist(timings))
  )

  saveRDS(results, output_file)
  invisible(results)
}

print_timing_comparison <- function(cpp, rust) {
  steps <- names(cpp$timings)
  cat("\n==> Timing comparison\n")
  cat(sprintf("%-24s %10s %10s %12s\n", "Step", "Seurat", "RSeurat", "Speedup"))
  cat(strrep("-", 62), "\n", sep = "")
  for (step in steps) {
    cpp_time <- cpp$timings[[step]]
    rust_time <- rust$timings[[step]]
    speedup <- if (rust_time > 0) cpp_time / rust_time else NA_real_
    cat(sprintf("%-24s %10.3f %10.3f %11.2fx\n", step, cpp_time, rust_time, speedup))
  }
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf(
    "%-24s %10.3f %10.3f %11.2fx\n",
    "Total native kernels",
    cpp$native_kernel_seconds,
    rust$native_kernel_seconds,
    cpp$native_kernel_seconds / rust$native_kernel_seconds
  ))
  cat(sprintf(
    "%-24s %10.3f %10.3f %11.2fx\n",
    "Total workflow",
    cpp$total_seconds,
    rust$total_seconds,
    cpp$total_seconds / rust$total_seconds
  ))
}

compare_outputs <- function(cpp, rust) {
  exact_fields <- c(
    "n_cells",
    "n_genes",
    "n_clusters",
    "cluster_table",
    "variable_features_digest",
    "normalized_digest",
    "scaled_rounded_digest",
    "snn_digest",
    "snn_nnz",
    "cluster_digest"
  )
  informational_fields <- c("scaled_digest", "umap_digest")

  cat("\n==> Output parity\n")
  all_ok <- TRUE
  for (field in exact_fields) {
    ok <- identical(cpp[[field]], rust[[field]])
    all_ok <- all_ok && ok
    cat(sprintf("  %-26s %s\n", field, if (ok) "OK" else "MISMATCH"))
    if (!ok) {
      cat("    Seurat :", paste(cpp[[field]], collapse = ", "), "\n")
      cat("    RSeurat:", paste(rust[[field]], collapse = ", "), "\n")
    }
  }
  for (field in informational_fields) {
    ok <- identical(cpp[[field]], rust[[field]])
    cat(sprintf(
      "  %-26s %s\n",
      field,
      if (ok) "OK" else "MISMATCH (informational)"
    ))
  }

  if (!all_ok) {
    stop("Output mismatch between Seurat C++ and RSeurat Rust runs.", call. = FALSE)
  }
}

cat("==> Upstream workflow: ", UPSTREAM_WORKFLOW_URL, "\n", sep = "")
cat("==> Data source: ", PBMC3K_DATA_URL, "\n\n", sep = "")

download_pbmc3k()
counts <- Read10X(data.dir = pbmc_matrix_dir)
# Seurat replaces underscores in feature names internally; normalize up front to avoid warnings.
rownames(counts) <- gsub("_", "-", rownames(counts))

cpp_file <- file.path(output_root, "seurat_cpp_results.rds")
rust_file <- file.path(output_root, "rseurat_rust_results.rds")

cat("==> Running seurat-standard-analysis PBMC workflow with Seurat C++ kernels...\n")
cpp <- run_workflow(make_backend("cpp"), counts, cpp_file)

cat("\n==> Running seurat-standard-analysis PBMC workflow with RSeurat Rust kernels...\n")
rust <- run_workflow(make_backend("rust"), counts, rust_file)

compare_outputs(cpp, rust)
print_timing_comparison(cpp, rust)

cat("\nAll exact parity checks passed. Speedup > 1.0 means RSeurat Rust was faster.\n")
cat("Results saved in: ", output_root, "\n", sep = "")
