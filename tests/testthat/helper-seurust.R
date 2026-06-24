#' Skip when seurust is not installed (Rust backend package).
#' @keywords internal
skip_if_no_seurust <- function() {
  if (!requireNamespace("seurust", quietly = TRUE)) {
    skip("seurust not installed; install the sibling package from seurust/")
  }
}
