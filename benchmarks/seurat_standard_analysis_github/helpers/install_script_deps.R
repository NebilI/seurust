# Install R packages required by ajtimon/seurat-standard-analysis scripts.
install_upstream_script_deps <- function() {
  repos <- c(CRAN = "https://cloud.r-project.org")
  options(repos = repos)

  cran_pkgs <- c(
    "pacman", "dplyr", "patchwork", "scales", "future", "future.apply",
    "parallelDist", "viridis", "pheatmap", "janitor", "purrr", "tidyr",
    "igraph", "ggalluvial"
  )
  for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("==> Installing ", pkg, "...")
      install.packages(pkg, quiet = TRUE)
    }
  }

  if (!requireNamespace("harmony", quietly = TRUE)) {
    message("==> Installing harmony...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = repos, quiet = TRUE)
    }
    tryCatch(
      install.packages("harmony", repos = repos, quiet = TRUE),
      error = function(e) NULL
    )
    if (!requireNamespace("harmony", quietly = TRUE)) {
      remotes::install_github(
        "immunogenomics/harmony",
        upgrade = "never",
        quiet = TRUE
      )
    }
  }
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Failed to install harmony for 02_gbm_seurat_adapted.R.", call. = FALSE)
  }

  if (!requireNamespace("presto", quietly = TRUE)) {
    message("==> Installing presto from GitHub...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", quiet = TRUE)
    }
    remotes::install_github("immunogenomics/presto", upgrade = "never", quiet = TRUE)
  }

  bioc_pkgs <- c("SingleR", "celldex", "msigdbr")
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_bioc) > 0) {
    message("==> Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", quiet = TRUE)
    }
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
  }

  optional_pkgs <- c("gprofiler2", "geomtextpath", "ggnewscale", "readxl", "cluster", "mclust")
  for (pkg in optional_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch(
        install.packages(pkg, quiet = TRUE),
        error = function(e) {
          message("Optional package not installed: ", pkg)
        }
      )
    }
  }

  invisible(TRUE)
}
