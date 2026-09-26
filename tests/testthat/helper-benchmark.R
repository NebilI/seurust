#' Time a pair of Seurat (C++) and seurust callables with warmup and repeated runs.
#'
#' @param cpp_fn Zero-argument function calling Seurat's C++ backend (`Seurat:::`).
#' @param rust_fn Zero-argument function calling seurust (`seurust::`).
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
  rust <- seurust::ComputeSNN(nn, prune)
  if (!isTRUE(all.equal(as.matrix(rust), as.matrix(cpp), tolerance = 1e-10))) {
    stop("ComputeSNN parity failed for ", label, call. = FALSE)
  }
  bench <- benchmark_rust_cpp(
    cpp_fn = function() Seurat:::ComputeSNN(nn, prune),
    rust_fn = function() seurust::ComputeSNN(nn, prune),
    n_warmup = n_warmup,
    n_reps = n_reps
  )
  attr(bench, "label") <- label
  bench
}

#' Resident set size of this R process in megabytes (Linux).
#' @keywords internal
process_rss_mb <- function() {
  if (!file.exists("/proc/self/statm")) {
    return(NA_real_)
  }
  fields <- scan("/proc/self/statm", quiet = TRUE, what = numeric(), nmax = 2)
  if (length(fields) < 2L) {
    return(NA_real_)
  }
  page <- 4096
  conf <- suppressWarnings(
    system2("getconf", "PAGESIZE", stdout = TRUE, stderr = FALSE)
  )
  conf_n <- suppressWarnings(as.numeric(conf))
  if (length(conf_n) == 1L && is.finite(conf_n) && conf_n > 0) {
    page <- conf_n
  }
  unname(fields[[2]] * page / (1024^2))
}

#' R-heap megabytes from a `gc()` matrix (`used` or `max used`).
#' `gc()` labels cell counts as `used` / `max used` and megabytes as a following `Mb` column.
#' @keywords internal
gc_mb_column <- function(g, which = c("used", "max")) {
  which <- match.arg(which)
  cn <- colnames(g)
  mb_cols <- which(cn %in% c("Mb", "(Mb)"))
  idx <- if (identical(which, "used")) 1L else 3L
  if (length(mb_cols) >= idx) {
    return(unname(sum(g[, mb_cols[[idx]]])))
  }
  col <- if (identical(which, "used")) 2L else 6L
  if (ncol(g) >= col) {
    return(unname(sum(g[, col])))
  }
  NA_real_
}

gc_used_mb <- function(g = gc()) {
  gc_mb_column(g, "used")
}

gc_max_used_mb <- function(g = gc()) {
  gc_mb_column(g, "max")
}

#' Measure R-heap peak and process RSS while holding a function result.
#' @keywords internal
measure_memory_fn <- function(fn, n_reps = 5L) {
  n_reps <- as.integer(n_reps)
  stopifnot(n_reps >= 1L)
  rss <- rep(NA_real_, n_reps)
  heap <- rep(NA_real_, n_reps)
  result_bytes <- rep(NA_real_, n_reps)
  for (i in seq_len(n_reps)) {
    g0 <- gc(reset = TRUE, full = TRUE)
    used0 <- gc_used_mb(g0)
    rss0 <- process_rss_mb()
    out <- fn()
    rss1 <- process_rss_mb()
    g1 <- gc()
    rss[i] <- if (is.na(rss0) || is.na(rss1)) NA_real_ else max(rss1 - rss0, 0)
    heap[i] <- max(gc_max_used_mb(g1) - used0, 0)
    result_bytes[i] <- as.numeric(utils::object.size(out))
    rm(out)
  }
  list(
    rss_delta_mb = suppressWarnings(stats::median(rss, na.rm = TRUE)),
    r_heap_mb = suppressWarnings(stats::median(heap, na.rm = TRUE)),
    result_mb = suppressWarnings(stats::median(result_bytes, na.rm = TRUE)) / (1024^2)
  )
}

#' Time and measure memory for a Seurat (C++) vs seurust pair.
#' @keywords internal
compare_rust_cpp <- function(
    cpp_fn,
    rust_fn,
    n_warmup = 1L,
    n_reps = 10L,
    n_mem = 5L) {
  bench <- benchmark_rust_cpp(
    cpp_fn = cpp_fn,
    rust_fn = rust_fn,
    n_warmup = n_warmup,
    n_reps = n_reps
  )
  cpp_mem <- measure_memory_fn(cpp_fn, n_reps = n_mem)
  rust_mem <- measure_memory_fn(rust_fn, n_reps = n_mem)
  heap_ratio <- if (is.finite(rust_mem$r_heap_mb) && rust_mem$r_heap_mb > 0) {
    cpp_mem$r_heap_mb / rust_mem$r_heap_mb
  } else {
    NA_real_
  }
  rss_ratio <- if (is.finite(rust_mem$rss_delta_mb) && rust_mem$rss_delta_mb > 0) {
    cpp_mem$rss_delta_mb / rust_mem$rss_delta_mb
  } else {
    NA_real_
  }
  c(
    bench,
    list(
      cpp_mem = cpp_mem,
      rust_mem = rust_mem,
      rust_vs_cpp_mem = unname(heap_ratio),
      rust_vs_cpp_rss = unname(rss_ratio)
    )
  )
}

#' Format a combined runtime + memory comparison line.
#' @keywords internal
format_compare <- function(cmp, label) {
  time_line <- format_benchmark(cmp, label)
  sprintf(
    paste0(
      "%s; ",
      "C++ heap=%.2f MB rssΔ=%.2f MB result=%.3f MB; ",
      "Rust heap=%.2f MB rssΔ=%.2f MB result=%.3f MB; ",
      "C++/Rust heap=%.2fx rss=%.2fx"
    ),
    time_line,
    cmp$cpp_mem$r_heap_mb,
    cmp$cpp_mem$rss_delta_mb,
    cmp$cpp_mem$result_mb,
    cmp$rust_mem$r_heap_mb,
    cmp$rust_mem$rss_delta_mb,
    cmp$rust_mem$result_mb,
    if (is.finite(cmp$rust_vs_cpp_mem)) cmp$rust_vs_cpp_mem else NA_real_,
    if (is.finite(cmp$rust_vs_cpp_rss)) cmp$rust_vs_cpp_rss else NA_real_
  )
}

#' Print a runtime+memory comparison and assert measurements are finite.
#' @keywords internal
expect_compare_report <- function(cmp, label) {
  line <- format_compare(cmp, label)
  cat(line, "\n", sep = "")
  testthat::expect_true(
    is.finite(cmp$rust_vs_cpp) && cmp$rust_vs_cpp > 0,
    info = line
  )
  testthat::expect_true(
    is.finite(cmp$cpp_mem$r_heap_mb) && is.finite(cmp$rust_mem$r_heap_mb),
    info = paste0(line, " (R heap must be measurable)")
  )
  invisible(cmp)
}

.seurust_compare_env <- new.env(parent = emptyenv())
.seurust_compare_env$rows <- list()

#' Record one function's correctness + runtime + memory comparison.
#' @keywords internal
register_compare_row <- function(
    name,
    size,
    parity,
    cmp) {
  row <- list(
    function_name = as.character(name),
    size = as.character(size),
    parity = isTRUE(parity),
    n_reps = as.integer(cmp$n_reps),
    cpp_median_us = unname(cmp$cpp$median),
    rust_median_us = unname(cmp$rust$median),
    speedup_cpp_over_rust = unname(cmp$rust_vs_cpp),
    cpp_heap_mb = unname(cmp$cpp_mem$r_heap_mb),
    rust_heap_mb = unname(cmp$rust_mem$r_heap_mb),
    heap_ratio_cpp_over_rust = unname(cmp$rust_vs_cpp_mem),
    cpp_rss_delta_mb = unname(cmp$cpp_mem$rss_delta_mb),
    rust_rss_delta_mb = unname(cmp$rust_mem$rss_delta_mb),
    rss_ratio_cpp_over_rust = unname(cmp$rust_vs_cpp_rss),
    cpp_result_mb = unname(cmp$cpp_mem$result_mb),
    rust_result_mb = unname(cmp$rust_mem$result_mb)
  )
  .seurust_compare_env$rows[[length(.seurust_compare_env$rows) + 1L]] <- row
  invisible(row)
}

#' Collected comparison rows as a data frame.
#' @keywords internal
compare_results_df <- function() {
  rows <- .seurust_compare_env$rows
  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, lapply(rows, function(r) {
    as.data.frame(r, stringsAsFactors = FALSE)
  }))
}

reset_compare_results <- function() {
  .seurust_compare_env$rows <- list()
  invisible(NULL)
}

#' Default output paths for the merge-visible comparison report.
#' @keywords internal
compare_report_paths <- function() {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  md_env <- Sys.getenv("SEURUST_COMPARE_OUT", unset = "")
  csv_env <- Sys.getenv("SEURUST_COMPARE_CSV", unset = "")
  md <- if (nzchar(md_env)) {
    md_env
  } else if (nzchar(workspace)) {
    file.path(workspace, "seurust-compare.md")
  } else {
    file.path(tempdir(), "seurust-compare.md")
  }
  csv <- if (nzchar(csv_env)) {
    csv_env
  } else if (nzchar(workspace)) {
    file.path(workspace, "seurust-compare.csv")
  } else {
    file.path(tempdir(), "seurust-compare.csv")
  }
  list(md = md, csv = csv)
}

#' Markdown table of collected Seurat vs seurust comparisons.
#' @keywords internal
format_compare_markdown <- function(df = compare_results_df()) {
  if (!nrow(df)) {
    return("No seurust vs Seurat comparison rows were recorded.\n")
  }
  header <- paste(
    "| Function | Size | Parity | C++ median (µs) | Rust median (µs) | Speedup (C++/Rust) |",
    "C++ result (MB) | Rust result (MB) | C++ heap (MB) | Rust heap (MB) |",
    sep = " "
  )
  sep <- paste(rep("| ---", 10L), collapse = " ")
  sep <- paste0(sep, " |")
  rows <- vapply(
    X = seq_len(nrow(df)),
    FUN = function(i) {
      r <- df[i, ]
      sprintf(
        "| %s | %s | %s | %.2f | %.2f | %.2fx | %.4f | %.4f | %.2f | %.2f |",
        r$function_name,
        r$size,
        if (isTRUE(r$parity)) "pass" else "FAIL",
        r$cpp_median_us,
        r$rust_median_us,
        r$speedup_cpp_over_rust,
        r$cpp_result_mb,
        r$rust_result_mb,
        r$cpp_heap_mb,
        r$rust_heap_mb
      )
    },
    FUN.VALUE = character(1)
  )
  paste0(
    "Speedup is **C++ time ÷ Rust time**. Values **> 1.0 mean Rust is faster**. ",
    "Result (MB) is `object.size` of the return value. Heap (MB) is R's peak `gc()` usage during the call.\n\n",
    header, "\n", sep, "\n", paste(rows, collapse = "\n"), "\n"
  )
}

#' Write markdown + CSV comparison reports for CI / local inspection.
#' @keywords internal
write_compare_report <- function(df = compare_results_df()) {
  paths <- compare_report_paths()
  md <- paste0(
    "# seurust vs Seurat (runtime + memory)\n\n",
    format_compare_markdown(df)
  )
  dir.create(dirname(paths$md), recursive = TRUE, showWarnings = FALSE)
  writeLines(md, paths$md)
  if (nrow(df)) {
    utils::write.csv(df, paths$csv, row.names = FALSE)
  }
  cat("Wrote comparison report to ", paths$md, "\n", sep = "")
  if (file.exists(paths$csv)) {
    cat("Wrote comparison CSV to ", paths$csv, "\n", sep = "")
  }
  invisible(paths)
}
