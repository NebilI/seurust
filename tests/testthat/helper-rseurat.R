#' Skip when SeuratRust is not installed (Rust backend package).
#' @keywords internal
skip_if_no_seuratrust <- function() {
  if (!requireNamespace("SeuratRust", quietly = TRUE)) {
    skip("SeuratRust not installed; install the sibling package from SeuratRust/")
  }
}
