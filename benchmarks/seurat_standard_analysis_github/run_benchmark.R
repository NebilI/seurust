#!/usr/bin/env Rscript
# Time upstream ajtimon/seurat-standard-analysis scripts with Seurat C++ vs RSeurat Rust.
# Upstream scripts are sourced unchanged; RSeurat is activated by patching Seurat natives.

UPSTREAM_REPO <- "https://github.com/ajtimon/seurat-standard-analysis"
UPSTREAM_CLONE_URL <- "https://github.com/ajtimon/seurat-standard-analysis.git"

find_repo_root <- function() {
  candidates <- c(
    normalizePath(".", winslash = "/", mustWork = FALSE),
    normalizePath("..", winslash = "/", mustWork = FALSE),
    normalizePath("../..", winslash = "/", mustWork = FALSE),
    normalizePath("../../..", winslash = "/", mustWork = FALSE)
  )
  for (path in candidates) {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        dir.exists(file.path(path, "RSeurat"))) {
      return(path)
    }
  }
  Sys.getenv("SEURAT_PKG_ROOT", unset = "/workspace")
}

ensure_upstream_repo <- function(upstream_root) {
  script_dir <- file.path(upstream_root, "code")
  if (dir.exists(script_dir)) {
    return(invisible(upstream_root))
  }
  dir.create(dirname(upstream_root), recursive = TRUE, showWarnings = FALSE)
  message("==> Cloning upstream repo into ", upstream_root, " ...")
  status <- system2(
    "git",
    c("clone", "--depth", "1", UPSTREAM_CLONE_URL, upstream_root),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L) || !dir.exists(script_dir)) {
    stop("Failed to clone upstream repo: ", UPSTREAM_REPO, call. = FALSE)
  }
  invisible(upstream_root)
}

empty_native_timings <- function() {
  new.env(parent = emptyenv())
}

record_native_timing <- function(store, fn_name, elapsed) {
  if (exists(fn_name, envir = store, inherits = FALSE)) {
    entry <- get(fn_name, envir = store, inherits = FALSE)
  } else {
    entry <- list(calls = 0L, seconds = 0)
  }
  entry$calls <- entry$calls + 1L
  entry$seconds <- entry$seconds + elapsed
  assign(fn_name, entry, envir = store)
}

native_timing_summary <- function(store) {
  names <- ls(store, all.names = TRUE)
  if (length(names) == 0L) {
    return(data.frame(function_name = character(), calls = integer(), seconds = numeric()))
  }
  rows <- lapply(names, function(fn_name) {
    entry <- get(fn_name, envir = store, inherits = FALSE)
    data.frame(
      function_name = fn_name,
      calls = entry$calls,
      seconds = entry$seconds,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$seconds, decreasing = TRUE), , drop = FALSE]
}

wrap_patched_native_timings <- function(patched_fns, store) {
  seurat_ns <- asNamespace("Seurat")
  for (fn_name in patched_fns) {
    if (!exists(fn_name, envir = seurat_ns, inherits = FALSE)) {
      next
    }
    target <- get(fn_name, envir = seurat_ns, inherits = FALSE)
    wrapper <- local({
      name <- fn_name
      fn <- target
      function(...) {
        call_args <- as.list(match.call(expand.dots = FALSE))[-1L]
        timing <- system.time(result <- do.call(fn, call_args, envir = parent.frame()))
        record_native_timing(store, name, unname(timing[["elapsed"]]))
        result
      }
    })
    formals(wrapper) <- formals(target)
    if (bindingIsLocked(fn_name, seurat_ns)) {
      unlockBinding(fn_name, seurat_ns)
    }
    assign(fn_name, wrapper, envir = seurat_ns)
    lockBinding(fn_name, seurat_ns)
  }
  invisible(store)
}

run_timed_script <- function(script_path, engine, code_dir) {
  patched <- patch_seurat_backend(engine)
  native_timing_store <- empty_native_timings()
  wrap_patched_native_timings(patched, native_timing_store)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(code_dir)

  script_name <- basename(script_path)
  if (grepl("02_gbm", script_name, fixed = TRUE)) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("harmony is required for 02_gbm_seurat_adapted.R but is not installed.", call. = FALSE)
    }
    suppressPackageStartupMessages(library(harmony))
    enable_harmony_compat()
  }

  gc(verbose = FALSE)
  error_msg <- NULL
  timing <- system.time({
    err <- tryCatch(
      {
        source(basename(script_path), local = FALSE, echo = FALSE, chdir = FALSE)
        NULL
      },
      error = function(e) {
        error_msg <<- conditionMessage(e)
      }
    )
  })

  list(
    script = basename(script_path),
    backend = backend_label(engine),
    engine = engine,
    seconds = unname(timing[["elapsed"]]),
    native_timings = native_timing_summary(native_timing_store),
    success = is.null(error_msg),
    error = error_msg
  )
}

print_timing_table <- function(results) {
  scripts <- unique(vapply(results, `[[`, "", "script"))
  cat("\n==> Timing comparison (upstream scripts, wall-clock seconds)\n")
  cat(sprintf("%-32s %12s %12s %10s\n", "Script", "Seurat C++", "RSeurat", "Speedup"))
  cat(strrep("-", 72), "\n", sep = "")

  for (script in scripts) {
    rows <- Filter(function(x) identical(x$script, script), results)
    cpp <- rows[[which(vapply(rows, function(x) x$engine == "cpp", logical(1)))[1]]]
    rust <- rows[[which(vapply(rows, function(x) x$engine == "rust", logical(1)))[1]]]

    cpp_time <- if (isTRUE(cpp$success)) cpp$seconds else NA_real_
    rust_time <- if (isTRUE(rust$success)) rust$seconds else NA_real_
    speedup <- if (!is.na(cpp_time) && !is.na(rust_time) && rust_time > 0) {
      cpp_time / rust_time
    } else {
      NA_real_
    }

    cat(sprintf(
      "%-32s %12.1f %12.1f %9.2fx\n",
      script,
      cpp_time,
      rust_time,
      speedup
    ))

    if (!isTRUE(cpp$success)) {
      cat("  Seurat C++ failed: ", cpp$error, "\n", sep = "")
    }
    if (!isTRUE(rust$success)) {
      cat("  RSeurat failed: ", rust$error, "\n", sep = "")
    }
  }

  ok <- Filter(function(x) isTRUE(x$success), results)
  if (length(ok) >= 2L) {
    cpp_total <- sum(vapply(Filter(function(x) x$engine == "cpp", ok), `[[`, 0, "seconds"))
    rust_total <- sum(vapply(Filter(function(x) x$engine == "rust", ok), `[[`, 0, "seconds"))
    cat(strrep("-", 72), "\n", sep = "")
    cat(sprintf(
      "%-32s %12.1f %12.1f %9.2fx\n",
      "Total (successful scripts)",
      cpp_total,
      rust_total,
      if (rust_total > 0) cpp_total / rust_total else NA_real_
    ))
  }
}

native_seconds <- function(result) {
  timings <- result$native_timings
  if (is.null(timings) || nrow(timings) == 0L) {
    return(0)
  }
  sum(timings$seconds)
}

print_native_timing_table <- function(results) {
  scripts <- unique(vapply(results, `[[`, "", "script"))
  cat("\n==> Patched native timing breakdown (seconds)\n")
  cat(sprintf("%-32s %-12s %12s %12s %12s\n", "Script", "Backend", "Total", "Native", "Residual"))
  cat(strrep("-", 86), "\n", sep = "")
  for (script in scripts) {
    rows <- Filter(function(x) identical(x$script, script) && isTRUE(x$success), results)
    for (row in rows) {
      native <- native_seconds(row)
      residual <- max(0, row$seconds - native)
      cat(sprintf(
        "%-32s %-12s %12.1f %12.1f %12.1f\n",
        script,
        row$backend,
        row$seconds,
        native,
        residual
      ))
    }
  }
}

repo_root <- find_repo_root()
setwd(repo_root)

system2("Rscript", "docker/scripts/bootstrap-dev-env.R", stdout = "", stderr = "")

benchmark_root <- file.path(repo_root, "benchmarks/seurat_standard_analysis_github")
upstream_root <- Sys.getenv(
  "SEURAT_STANDARD_ANALYSIS_ROOT",
  unset = file.path(benchmark_root, "upstream")
)
output_root <- Sys.getenv(
  "SEURAT_STANDARD_ANALYSIS_OUTPUT_DIR",
  unset = file.path(benchmark_root, "output")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

helper_dir <- file.path(benchmark_root, "helpers")
source(file.path(helper_dir, "backend_patch.R"), local = TRUE)
source(file.path(helper_dir, "harmony_compat.R"), local = TRUE)
source(file.path(helper_dir, "output_parity.R"), local = TRUE)
source(file.path(helper_dir, "install_script_deps.R"), local = TRUE)

suppressPackageStartupMessages({
  devtools::load_all(recompile = FALSE, quiet = TRUE)
  library(Seurat)
  library(RSeurat)
})

ensure_upstream_repo(upstream_root)
install_upstream_script_deps()
enable_harmony_compat()

scripts <- c(
  "01_pbmc_satija_tutorial.R",
  "02_gbm_seurat_adapted.R",
  "03_gbmap_exploration.R"
)
code_dir <- file.path(upstream_root, "code")
missing <- scripts[!file.exists(file.path(code_dir, scripts))]
if (length(missing) > 0) {
  stop("Missing upstream scripts: ", paste(missing, collapse = ", "), call. = FALSE)
}

cat("==> Upstream repo: ", UPSTREAM_REPO, "\n", sep = "")
cat("==> Running upstream scripts unchanged; RSeurat enabled via Seurat native patch.\n\n")

results <- list()
parity_failures <- character()
for (script in scripts) {
  script_path <- file.path(code_dir, script)
  cpp_digests <- NULL
  for (engine in c("cpp", "rust")) {
    cat("==> ", script, " [", backend_label(engine), "]\n", sep = "")
    res <- run_timed_script(script_path, engine, code_dir)
    results[[length(results) + 1L]] <- res
    if (isTRUE(res$success)) {
      cat("    Finished in ", sprintf("%.1f", res$seconds), " s\n", sep = "")
      run_digests <- capture_script_outputs(script, code_dir)
      if (identical(engine, "cpp")) {
        cpp_digests <- run_digests
      } else if (!is.null(cpp_digests)) {
        parity_ok <- tryCatch(
          {
            compare_script_outputs(script, cpp_digests, run_digests)
            TRUE
          },
          error = function(e) {
            parity_failures <<- c(parity_failures, script)
            cat("  ", conditionMessage(e), "\n", sep = "")
            FALSE
          }
        )
        if (!parity_ok) {
          res$parity_ok <- FALSE
          results[[length(results)]] <- res
        }
      }
    } else {
      cat("    FAILED after ", sprintf("%.1f", res$seconds), " s: ", res$error, "\n", sep = "")
    }
  }
  cat("\n")
}

saveRDS(results, file.path(output_root, "script_timing_results.rds"))
print_timing_table(results)
print_native_timing_table(results)

if (length(parity_failures) > 0L) {
  stop(
    "Output parity check failed for: ",
    paste(unique(parity_failures), collapse = ", "),
    call. = FALSE
  )
}

cat("\nAll successful script pairs matched Seurat C++ and RSeurat outputs.\n")
cat("Speedup > 1.0 means RSeurat was faster for that script.\n")
cat("Detailed results: ", file.path(output_root, "script_timing_results.rds"), "\n", sep = "")
