#!/usr/bin/env bash
# Build Seurat (C++) and SeuratRust, then run cross-package parity checks.
set -euo pipefail

cd /workspace

export NOT_CRAN="${NOT_CRAN:-1}"
export SEURAT_KEEP_RUST_TARGET="${SEURAT_KEEP_RUST_TARGET:-1}"

echo "==> Installing Seurat Depends/Imports (needed for compile_dll)..."
Rscript docker/scripts/install-imports.R

echo "==> Installing Seurat (C++/Rcpp only)..."
Rscript -e "pkgbuild::compile_dll(debug = FALSE, compile_attributes = FALSE)"

echo "==> Installing SeuratRust..."
# Windows checkouts may have CRLF in shell/config files; strip before R CMD INSTALL.
sed -i 's/\r$//' SeuratRust/configure SeuratRust/cleanup SeuratRust/DESCRIPTION SeuratRust/src/entrypoint.c 2>/dev/null || true
export NOT_CRAN=1 SEURAT_KEEP_RUST_TARGET=1
cd SeuratRust
Rscript tools/config.R
cd ..
R CMD INSTALL --preclean SeuratRust

echo "==> Running parity checks..."
Rscript docker/scripts/run-rust-parity.R

echo "==> Done."
