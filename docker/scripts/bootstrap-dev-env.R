#!/usr/bin/env Rscript
# Ensure Seurat dev DLL, RSeurat, and Imports exist (each `docker compose run`
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

rseurat_ok <- function() {
  if (!requireNamespace("RSeurat", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch(
    {
      x <- c(1, 2)
      i <- c(0L, 1L)
      out <- RSeurat::row_sum_dgcmatrix(x, i, 2L, 2L)
      length(out) == 2L
    },
    error = function(e) FALSE
  )
}

if (!rseurat_ok()) {
  message("==> Installing RSeurat...")
  run_or_stop(
    "bash",
    c(
      "-c",
      paste(
        "export NOT_CRAN=1 SEURAT_KEEP_RUST_TARGET=1;",
        "sed -i 's/\\r$//' RSeurat/configure RSeurat/cleanup",
        "RSeurat/DESCRIPTION RSeurat/src/entrypoint.c 2>/dev/null || true;",
        "cd RSeurat && Rscript tools/config.R && cd ..;",
        "R CMD INSTALL --preclean RSeurat"
      )
    )
  )
}

message("==> Dev environment ready.")
