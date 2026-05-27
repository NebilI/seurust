#!/usr/bin/env bash
# Build the Seurat DLL (Rcpp + Rust) and compare C++ vs Rust implementations.
set -euo pipefail

cd /workspace

export NOT_CRAN="${NOT_CRAN:-1}"
export SEURAT_KEEP_RUST_TARGET="${SEURAT_KEEP_RUST_TARGET:-1}"

echo "==> Configuring Makevars (Rust + Rcpp)..."
Rscript tools/config.R

echo "==> Regenerating Rcpp exports (unified DLL init)..."
Rscript tools/fix-rcpp-init.R

echo "==> Compiling shared library..."
Rscript -e "pkgbuild::compile_dll(debug = FALSE, compile_attributes = FALSE)"

echo "==> Running parity checks and testthat files..."
Rscript docker/scripts/run-rust-parity.R

echo "==> Done."
