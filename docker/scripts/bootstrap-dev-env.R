#!/usr/bin/env Rscript
# Ensure seurust is installed inside a fresh container. Seurat is installed
# from CRAN when it is missing so parity examples can call Seurat:::.

setwd("/workspace")

run_or_stop <- function(cmd, args = character()) {
  status <- system2(cmd, args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Command failed: ", cmd, " ", paste(args, collapse = " "))
  }
  invisible(status)
}

if (!requireNamespace("Matrix", quietly = TRUE) || !requireNamespace("Seurat", quietly = TRUE)) {
  message("==> Installing Matrix and Seurat from CRAN...")
  install.packages(
    c("Matrix", "Seurat"),
    repos = c("https://cloud.r-project.org", "https://satijalab.r-universe.dev")
  )
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
