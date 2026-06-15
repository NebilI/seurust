#' Time a pair of Seurat (C++) and RSeurat callables with warmup and repeated runs.
#'
#' @param cpp_fn Zero-argument function calling Seurat's C++ backend (`Seurat:::`).
#' @param rust_fn Zero-argument function calling RSeurat (`RSeurat::`).
#' @param n_warmup Warmup iterations (not timed).
#' @param n_reps Timed repetitions; mean, sd, and median elapsed microseconds are reported.
#' @return A list with per-backend summaries, raw `times` vectors (microseconds), and
#'   `rust_vs_cpp` ratio from medians (>1 means Rust is faster).
#' @keywords internal
benchmark_rust_cpp <- function(cpp_fn, rust_fn, n_warmup = 3L, n_reps = 100L) {
  n_warmup <- as.integer(n_warmup)
  n_reps <- as.integer(n_reps)
  stopifnot(n_warmup >= 0L, n_reps >= 1L)

  for (w in seq_len(n_warmup)) {
    invisible(cpp_fn())
    invisible(rust_fn())
  }

  time_fn <- function(fn) {
    times_us <- if (requireNamespace("microbenchmark", quietly = TRUE)) {
      as.numeric(
        microbenchmark::microbenchmark(
          fn(),
          times = n_reps,
          warmup = 0L
        )$time
      ) / 1000
    } else {
      vapply(
        X = seq_len(n_reps),
        FUN = function(i) {
          t0 <- proc.time()[["elapsed"]]
          fn()
          (proc.time()[["elapsed"]] - t0) * 1e6
        },
        FUN.VALUE = numeric(1)
      )
    }
    list(
      times = times_us,
      median = stats::median(times_us),
      mean = mean(times_us),
      sd = stats::sd(times_us),
      min = min(times_us),
      max = max(times_us)
    )
  }

  cpp <- time_fn(cpp_fn)
  rust <- time_fn(rust_fn)
  # Medians can be 0 when timer resolution exceeds runtime; fall back to means.
  cpp_basis <- if (cpp$median > 0) cpp$median else cpp$mean
  rust_basis <- if (rust$median > 0) rust$median else rust$mean
  if (rust_basis <= 0) {
    rust_basis <- max(rust$min, .Machine$double.eps)
  }
  if (cpp_basis <= 0) {
    cpp_basis <- max(cpp$min, .Machine$double.eps)
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
      "C++ mean=%.2f us (sd=%.2f), median=%.2f us; ",
      "Rust mean=%.2f us (sd=%.2f), median=%.2f us; ",
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
expect_rust_faster <- function(bench, label, tolerance = 0.95) {
  msg <- format_benchmark(bench, label)
  testthat::expect_true(
    bench$rust_vs_cpp >= tolerance,
    info = paste0(msg, " (goal: Rust/C++ median ratio >= ", tolerance, ")")
  )
}

#' Synthetic ranked-neighbor matrix for ComputeSNN benchmarks.
#' @keywords internal
make_compute_snn_nn <- function(n_cells, k = 20L, seed = 1L) {
  set.seed(seed)
  nn <- matrix(
    sample.int(n_cells, n_cells * k, replace = TRUE),
    nrow = n_cells,
    ncol = k
  )
  storage.mode(nn) <- "double"
  nn
}

#' Parity-check and time ComputeSNN for a given cell count.
#' @keywords internal
benchmark_compute_snn <- function(
    n_cells,
    k = 20L,
    prune = 0.01,
    label = NULL,
    n_warmup = 5L,
    n_reps = 100L,
    seed = 1L) {
  if (is.null(label)) {
    label <- sprintf("ComputeSNN (%d cells, k=%d)", n_cells, k)
  }
  nn <- make_compute_snn_nn(n_cells = n_cells, k = k, seed = seed)
  cpp <- Seurat:::ComputeSNN(nn, prune)
  rust <- RSeurat::ComputeSNN(nn, prune)
  if (!isTRUE(all.equal(as.matrix(rust), as.matrix(cpp), tolerance = 1e-10))) {
    stop("ComputeSNN parity failed for ", label, call. = FALSE)
  }
  bench <- benchmark_rust_cpp(
    cpp_fn = function() Seurat:::ComputeSNN(nn, prune),
    rust_fn = function() RSeurat::ComputeSNN(nn, prune),
    n_warmup = n_warmup,
    n_reps = n_reps
  )
  attr(bench, "label") <- label
  bench
}
