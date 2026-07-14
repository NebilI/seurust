source("tools/msrv.R")

env_debug <- Sys.getenv("DEBUG")
env_not_cran <- Sys.getenv("NOT_CRAN")

vendor_exists <- file.exists("src/rust/vendor.tar.xz")

is_not_cran <- env_not_cran != ""
is_debug <- env_debug != ""

if (is_debug) {
  is_not_cran <- TRUE
  message("Creating DEBUG build.")
}

if (!is_not_cran) {
  message("Building for CRAN.")
}

.cran_flags <- ifelse(
  !is_not_cran && vendor_exists,
  "-j 2 --offline",
  ""
)

.profile <- ifelse(is_debug, "", "--release")
.keep_rust_target <- {
  val <- Sys.getenv("SEURAT_KEEP_RUST_TARGET", unset = "")
  if (nzchar(val)) {
    !identical(tolower(val), "0") && !identical(tolower(val), "false")
  } else {
    is_not_cran || is_debug
  }
}
.clean_targets <- ifelse(is_debug || .keep_rust_target, "", "$(TARGET_DIR)")

webr_target <- "wasm32-unknown-emscripten"
is_wasm <- identical(R.version$platform, webr_target)

if (is_wasm) {
  message("Building for WebR")
}

target_libpath <- if (is_wasm) "wasm32-unknown-emscripten" else NULL
cfg <- if (is_debug) "debug" else "release"

.libdir <- paste(c(target_libpath, cfg), collapse = "/")

.target <- ifelse(is_wasm, paste0("--target=", webr_target), "")

.panic_exports <- ifelse(
  is_wasm,
  "CARGO_PROFILE_DEV_PANIC=\"abort\" CARGO_PROFILE_RELEASE_PANIC=\"abort\" ",
  ""
)

# Dev builds may regenerate R/extendr-wrappers.R via the document binary.
# CRAN/source installs ship committed wrappers and only build --lib.
.run_document <- if (is_not_cran) "true" else "false"

is_windows <- .Platform[["OS.type"]] == "windows"

mv_fp <- ifelse(
  is_windows,
  "src/Makevars.win.in",
  "src/Makevars.in"
)

mv_ofp <- ifelse(
  is_windows,
  "src/Makevars.win",
  "src/Makevars"
)

if (file.exists(mv_ofp)) {
  message("Cleaning previous `", mv_ofp, "`.")
  invisible(file.remove(mv_ofp))
}

mv_txt <- readLines(mv_fp)

new_txt <- gsub("@CRAN_FLAGS@", .cran_flags, mv_txt) |>
  gsub("@PROFILE@", .profile, x = _) |>
  gsub("@CLEAN_TARGET@", .clean_targets, x = _) |>
  gsub("@LIBDIR@", .libdir, x = _) |>
  gsub("@TARGET@", .target, x = _) |>
  gsub("@PANIC_EXPORTS@", .panic_exports, x = _) |>
  gsub("@RUN_DOCUMENT@", .run_document, x = _)

message("Writing `", mv_ofp, "`.")
con <- file(mv_ofp, open = "wb")
writeLines(new_txt, con, sep = "\n")
close(con)

message("`tools/config.R` has finished.")
