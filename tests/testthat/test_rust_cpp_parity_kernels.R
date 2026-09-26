# Parity tests for every ported seurust kernel versus Seurat's original C++ API.
# Complements the SNN/integration and modularity files with the remaining natives.

context("seurust/Seurat parity: all native kernels")

test_that("every Seurat native export has a seurust counterpart in the comparison catalog", {
  skip_if_no_seurust()
  cases <- seurust_kernel_cases()
  catalog_names <- vapply(cases, function(x) x$name, character(1))
  expect_equal(sort(catalog_names), sort(SEURAT_NATIVE_FUNCTIONS))
  seurust_exports <- getNamespaceExports("seurust")
  expect_true(all(SEURAT_NATIVE_FUNCTIONS %in% seurust_exports))
  seurat_ns <- getNamespace("Seurat")
  expect_true(all(vapply(
    SEURAT_NATIVE_FUNCTIONS,
    function(nm) exists(nm, envir = seurat_ns, inherits = FALSE),
    logical(1)
  )))
})

local({
  skip_if_no_seurust()
  cases <- seurust_kernel_cases()
  for (case in cases) {
    local({
      current <- case
      test_that(paste("seurust", current$name, "matches Seurat C++"), {
        skip_if_no_seurust()
        run_kernel_parity(current)
      })
    })
  }
})
