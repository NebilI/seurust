#!/usr/bin/env Rscript
# Ensure Seurat dev DLL, seurust, and Imports exist (each `docker compose run`
# starts a fresh container without prior install state).

pkg_root <- if (file.exists("DESCRIPTION")) {
  normalizePath(".", winslash = "/")
} else if (file.exists("../DESCRIPTION")) {
  normalizePath("..", winslash = "/")
} else {
  Sys.getenv("SEURAT_PKG_ROOT", unset = "/workspace")
}
setwd(pkg_root)

run_or_stop <- function(cmd, args = character()) {
  status <- system2(cmd, args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Command failed: ", cmd, " ", paste(args, collapse = " "))
  }
  invisible(status)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("==> Installing Seurat Imports...")
  run_or_stop("Rscript", "docker/scripts/install-imports.R")
}

if (length(Sys.glob("src/Seurat.{so,dll}")) == 0) {
  message("==> Compiling Seurat (C++)...")
  suppressPackageStartupMessages({
    library(pkgbuild)
  })
  compile_dll(debug = FALSE, compile_attributes = FALSE)
}

seurust_ok <- function() {
  if (!requireNamespace("seurust", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch(
    {
      x <- c(1, 2)
      i <- c(0L, 1L)
      out <- seurust::row_sum_dgcmatrix(x, i, 2L, 2L)
      length(out) == 2L
    },
    error = function(e) FALSE
  )
}

if (!seurust_ok()) {
  message("==> Installing seurust...")
  run_or_stop(
    "bash",
    c(
      "-c",
      paste(
        "export NOT_CRAN=1 SEURAT_KEEP_RUST_TARGET=1;",
        "sed -i 's/\\r$//' seurust/configure seurust/cleanup",
        "seurust/DESCRIPTION seurust/src/entrypoint.c 2>/dev/null || true;",
        "cd seurust && Rscript tools/config.R && cd ..;",
        "env NOT_CRAN=1 SEURAT_KEEP_RUST_TARGET=1 R CMD INSTALL --preclean seurust"
      )
    )
  )
}

message("==> Dev environment ready.")
