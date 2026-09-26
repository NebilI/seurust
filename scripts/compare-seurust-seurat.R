#!/usr/bin/env Rscript
# Compare every seurust kernel against original Seurat C++ for correctness,
# runtime, and memory. Writes seurust-compare.md and seurust-compare.csv.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  file.path(getwd(), "scripts/compare-seurust-seurat.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."))

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat must be installed to compare against the original library.")
  }
  if (!requireNamespace("seurust", quietly = TRUE)) {
    stop("seurust must be installed. From the repo root: R CMD INSTALL seurust")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix must be installed.")
  }
  library(Matrix)
  library(Seurat)
  library(seurust)
})

source(file.path(root, "tests/testthat/helper-seurust.R"), local = TRUE)
source(file.path(root, "tests/testthat/helper-benchmark.R"), local = TRUE)
source(file.path(root, "tests/testthat/helper-compare.R"), local = TRUE)

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required to run the comparison checks.")
}

reset_compare_results()
cases <- seurust_kernel_cases()
n_reps <- compare_n_reps()
n_mem <- compare_n_mem()
failed <- character()

cat("==> seurust vs Seurat: correctness, runtime, and memory\n")
cat("    ", length(cases), " kernels, ", n_reps, " timed reps, ", n_mem, " memory reps\n\n", sep = "")

for (case in cases) {
  out <- run_kernel_compare(case, n_warmup = 1L, n_reps = n_reps, n_mem = n_mem)
  label <- sprintf("%s [%s]", case$name, case$size)
  cat(format_compare(out$cmp, label), "\n", sep = "")
  if (!isTRUE(out$parity_ok)) {
    failed <- c(failed, paste0(case$name, ": ", out$parity_error))
  }
}

paths <- write_compare_report()
cat("\n", format_compare_markdown(), sep = "")

if (length(failed)) {
  stop(
    "Parity failed for:\n",
    paste0("  - ", failed, collapse = "\n"),
    call. = FALSE
  )
}

cat("\nAll ", length(cases), " kernels matched Seurat. Reports:\n", sep = "")
cat("  ", paths$md, "\n", sep = "")
if (file.exists(paths$csv)) {
  cat("  ", paths$csv, "\n", sep = "")
}
