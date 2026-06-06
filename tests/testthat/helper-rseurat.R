#' Skip when RSeurat is not installed (Rust backend package).
#' @keywords internal
skip_if_no_rseurat <- function() {
  if (!requireNamespace("RSeurat", quietly = TRUE)) {
    skip("RSeurat not installed; install the sibling package from RSeurat/")
  }
}
