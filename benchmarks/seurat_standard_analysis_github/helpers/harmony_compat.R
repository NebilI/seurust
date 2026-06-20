# Map deprecated harmony Seurat arguments to their current names.
# Upstream 02_gbm_seurat_adapted.R passes project.pca; current harmony expects project.dim.
# The upstream script also reinstalls/reloads harmony via pacman, so we re-apply on package load.

.harmony_compat_store <- new.env(parent = emptyenv())

.translate_harmony_args <- function(dots) {
  if ("project.pca" %in% names(dots)) {
    if (!"project.dim" %in% names(dots)) {
      dots$project.dim <- dots$project.pca
    }
    dots$project.pca <- NULL
  }
  dots
}

.current_harmony_namespace <- function() {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony is required for 02_gbm_seurat_adapted.R but is not installed.", call. = FALSE)
  }
  asNamespace("harmony")
}

.run_harmony_impl <- function(object, dots) {
  ns <- .current_harmony_namespace()
  if (inherits(object, "Seurat") && exists("RunHarmony.Seurat", envir = ns, inherits = FALSE)) {
    method <- get("RunHarmony.Seurat", envir = ns, inherits = FALSE)
    return(do.call(method, c(list(object = object), dots)))
  }
  generic <- get("RunHarmony", envir = ns, inherits = FALSE)
  do.call(generic, c(list(object), dots))
}

.run_harmony_compat <- function(object, ...) {
  dots <- .translate_harmony_args(list(...))
  .run_harmony_impl(object, dots)
}

.apply_global_harmony_compat_patch <- function() {
  assign("RunHarmony", .run_harmony_compat, envir = .GlobalEnv)
  invisible(TRUE)
}

.apply_harmony_compat_patch <- function() {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  ns <- asNamespace("harmony")
  if (!exists("RunHarmony", envir = ns, inherits = FALSE)) {
    stop("harmony does not export RunHarmony.", call. = FALSE)
  }

  # Do not patch the namespace generic: RunHarmony.Seurat calls the generic
  # internally to dispatch to RunHarmony.default, and wrapping it there breaks
  # that second dispatch. Patch only external lookup environments.
  .harmony_compat_store$original_fn <- get("RunHarmony", envir = ns, inherits = FALSE)

  package_env_name <- "package:harmony"
  if (package_env_name %in% search()) {
    package_env <- as.environment(package_env_name)
    if (bindingIsLocked("RunHarmony", package_env)) {
      unlockBinding("RunHarmony", package_env)
    }
    assign("RunHarmony", .run_harmony_compat, envir = package_env)
    lockBinding("RunHarmony", package_env)
  }

  invisible(TRUE)
}

.register_harmony_compat_hook <- function() {
  if (isTRUE(.harmony_compat_store$hook_registered)) {
    return(invisible(NULL))
  }
  setHook(packageEvent("harmony", "onLoad"), function(...) {
    .apply_harmony_compat_patch()
  })
  .harmony_compat_store$hook_registered <- TRUE
  invisible(NULL)
}

enable_harmony_compat <- function() {
  .register_harmony_compat_hook()
  .apply_harmony_compat_patch()
  .apply_global_harmony_compat_patch()
}
