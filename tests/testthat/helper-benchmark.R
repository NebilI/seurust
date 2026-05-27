#' Time a pair of C++ and Rust callables with warmup and repeated runs.
#'
#' @param cpp_fn Zero-argument function calling the C++ implementation.
#' @param rust_fn Zero-argument function calling the Rust implementation.
#' @param n_warmup Warmup iterations (not timed).
#' @param n_reps Timed repetitions; median elapsed seconds is reported.
#' @return A list with per-backend summaries and rust_vs_cpp ratio (>1 means Rust is faster).
#' @keywords internal
benchmark_rust_cpp <- function(cpp_fn, rust_fn, n_warmup = 3L, n_reps = 20L) {
  for (w in seq_len(n_warmup)) {
    invisible(cpp_fn())
    invisible(rust_fn())
  }

  time_fn <- function(fn) {
    times <- vapply(
      X = seq_len(n_reps),
      FUN = function(i) {
        t0 <- proc.time()[["elapsed"]]
        fn()
        proc.time()[["elapsed"]] - t0
      },
      FUN.VALUE = numeric(1)
    )
    c(
      median = stats::median(times),
      mean = mean(times),
      min = min(times),
      max = max(times)
    )
  }

  cpp <- time_fn(cpp_fn)
  rust <- time_fn(rust_fn)
  list(
    cpp = cpp,
    rust = rust,
    rust_vs_cpp = unname(cpp["median"] / rust["median"])
  )
}

#' Format benchmark output for logs / testthat messages.
#' @keywords internal
format_benchmark <- function(bench, label) {
  ratio <- bench$rust_vs_cpp
  winner <- if (ratio >= 1) {
    "Rust faster"
  } else {
    "C++ faster"
  }
  sprintf(
    "%s: C++ median=%.4fs, Rust median=%.4fs, Rust/C++=%.2fx (%s)",
    label,
    bench$cpp["median"],
    bench$rust["median"],
    ratio,
    winner
  )
}

#' Optionally fail when Rust is not faster than C++.
#' Set SEURAT_REQUIRE_RUST_FASTER=1 to enforce in CI or local runs.
#' @keywords internal
expect_rust_faster <- function(bench, label, tolerance = 1.0) {
  msg <- format_benchmark(bench, label)
  testthat::expect_true(
    bench$rust_vs_cpp >= tolerance,
    info = paste0(msg, " (goal: Rust >= ", tolerance, "x C++ speed)")
  )
}
