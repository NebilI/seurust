# Regenerate RcppExports.cpp with CallEntries + R_init, then rename the init
# function so entrypoint.c can provide the unified DLL initializer for Rcpp +
# extendr.
#
# pkgbuild::compile_dll() runs compileAttributes() while entrypoint.c is
# present, which omits R_init from RcppExports.cpp. Run this script before
# compile_dll(..., compile_attributes = FALSE).

pkg_root <- if (file.exists("DESCRIPTION")) {
  normalizePath(".", winslash = "/")
} else if (file.exists("../DESCRIPTION")) {
  normalizePath("..", winslash = "/")
} else {
  stop("Could not locate package root (missing DESCRIPTION).")
}

path <- file.path(pkg_root, "src/RcppExports.cpp")
entrypoint <- file.path(pkg_root, "src/entrypoint.c")
entry_bak <- file.path(pkg_root, "src/entrypoint.c.bak")

restore_entrypoint <- function() {
  if (file.exists(entry_bak)) {
    if (file.exists(entrypoint)) {
      file.remove(entrypoint)
    }
    file.rename(entry_bak, entrypoint)
  }
}

if (file.exists(entry_bak)) {
  restore_entrypoint()
}

if (file.exists(entrypoint)) {
  file.rename(entrypoint, entry_bak)
}

tryCatch(
  {
    setwd(pkg_root)
    Rcpp::compileAttributes()

    lines <- readLines(path)
    idx <- grep("^RcppExport void R_init_Seurat\\(", lines)
    if (length(idx) != 1L) {
      stop("Could not find R_init_Seurat definition in ", path)
    }

    lines[idx] <- sub("R_init_Seurat", "R_init_Seurat_rcpp", lines[idx])
    writeLines(lines, path)
    message("Prepared R_init_Seurat_rcpp in ", path)
  },
  finally = restore_entrypoint()
)

if (!file.exists(entrypoint)) {
  stop("Failed to restore ", entrypoint, " after regenerating Rcpp exports.")
}
