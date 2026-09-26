# Runtime and memory comparison for every ported seurust kernel vs Seurat C++.
# Results are printed and written to seurust-compare.md / .csv for merge CI.

context("seurust/Seurat runtime and memory")

test_that("all seurust kernels report runtime and memory versus Seurat", {
  skip_if_no_seurust()
  if (identical(Sys.getenv("SEURUST_SKIP_FULL_COMPARE"), "1")) {
    skip("Full runtime/memory table is published by the merge compare job")
  }
  reset_compare_results()
  cases <- seurust_kernel_cases()
  n_reps <- compare_n_reps()
  n_mem <- compare_n_mem()
  failures <- character()

  for (case in cases) {
    out <- run_kernel_compare(case, n_warmup = 1L, n_reps = n_reps, n_mem = n_mem)
    label <- sprintf("%s [%s]", case$name, case$size)
    cat(format_compare(out$cmp, label), "\n", sep = "")
    if (!isTRUE(out$parity_ok)) {
      failures <- c(failures, paste0(case$name, ": ", out$parity_error))
    }
    expect_true(
      is.finite(out$cmp$rust_vs_cpp) && out$cmp$rust_vs_cpp > 0,
      info = paste("non-finite speedup for", case$name)
    )
    expect_true(
      is.finite(out$cmp$cpp_mem$r_heap_mb) && is.finite(out$cmp$rust_mem$r_heap_mb),
      info = paste("R heap not measurable for", case$name)
    )
    if (identical(Sys.getenv("SEURAT_REQUIRE_RUST_FASTER"), "1")) {
      expect_rust_faster(out$cmp, label)
    }
  }

  paths <- write_compare_report()
  cat("\n", format_compare_markdown(), sep = "")
  expect_true(file.exists(paths$md))
  expect_length(failures, 0)
})
