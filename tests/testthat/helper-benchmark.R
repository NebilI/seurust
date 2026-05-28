#' Time a pair of Seurat (C++) and SeuratRust callables with warmup and repeated runs.
#'
#' @param cpp_fn Zero-argument function calling Seurat's C++ backend (`Seurat:::`).
#' @param rust_fn Zero-argument function calling SeuratRust (`SeuratRust::`).
#' @param n_warmup Warmup iterations (not timed).
#' @param n_reps Timed repetitions; mean, sd, and median elapsed seconds are reported.
#' @return A list with per-backend summaries, raw `times` vectors, and `rust_vs_cpp`
#'   ratio from medians (>1 means Rust is faster).
#' @keywords internal
benchmark_rust_cpp <- function(cpp_fn, rust_fn, n_warmup = 3L, n_reps = 20L) {
  n_warmup <- as.integer(n_warmup)
  n_reps <- as.integer(n_reps)
  stopifnot(n_warmup >= 0L, n_reps >= 1L)

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
    list(
      times = times,
      median = stats::median(times),
      mean = mean(times),
      sd = stats::sd(times),
      min = min(times),
      max = max(times)
    )
  }

  cpp <- time_fn(cpp_fn)
  rust <- time_fn(rust_fn)
  # Medians can be 0 when proc.time() resolution exceeds runtime; fall back to means.
  cpp_basis <- if (cpp$median > 0) cpp$median else cpp$mean
  rust_basis <- if (rust$median > 0) rust$median else rust$mean
  if (rust_basis <= 0) {
    rust_basis <- max(rust$min, .Machine$double.eps)
  }
  list(
    n_reps = n_reps,
    cpp = cpp,
    rust = rust,
    rust_vs_cpp = unname(cpp_basis / rust_basis)
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
    paste0(
      "%s (n=%d): ",
      "C++ mean=%.4fs (sd=%.4f), median=%.4fs; ",
      "Rust mean=%.4fs (sd=%.4f), median=%.4fs; ",
      "Rust/C++=%.2fx (%s)"
    ),
    label,
    bench$n_reps,
    bench$cpp$mean,
    bench$cpp$sd,
    bench$cpp$median,
    bench$rust$mean,
    bench$rust$sd,
    bench$rust$median,
    ratio,
    winner
  )
}

#' Run timing benchmark, print to stdout, and register a testthat expectation.
#' @keywords internal
expect_timing_report <- function(bench, label) {
  line <- format_benchmark(bench, label)
  cat(line, "\n", sep = "")
  testthat::expect_true(
    is.finite(bench$rust_vs_cpp) && bench$rust_vs_cpp > 0,
    info = line
  )
  invisible(bench)
}

#' Optionally fail when Rust is not faster than C++.
#' Set SEURAT_REQUIRE_RUST_FASTER=1 to enforce in CI or local runs.
#' @keywords internal
expect_rust_faster <- function(bench, label, tolerance = 1.0) {
  msg <- format_benchmark(bench, label)
  testthat::expect_true(
    bench$rust_vs_cpp >= tolerance,
    info = paste0(msg, " (goal: Rust/C++ median ratio >= ", tolerance, ")")
  )
}
