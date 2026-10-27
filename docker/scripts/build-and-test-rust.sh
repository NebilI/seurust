#!/usr/bin/env bash
# Configure, install, and test the seurust package.
set -euo pipefail

cd /workspace

export NOT_CRAN="${NOT_CRAN:-1}"
export SEURAT_KEEP_RUST_TARGET="${SEURAT_KEEP_RUST_TARGET:-1}"
export PATH="${PATH}:${HOME}/.cargo/bin:/usr/local/cargo/bin"

echo "==> Installing seurust Depends/Imports..."
Rscript docker/scripts/install-imports.R

if ! Rscript -e 'quit(status = if (requireNamespace("testthat", quietly = TRUE)) 0 else 1)'; then
  echo "==> Installing testthat..."
  Rscript -e 'install.packages("testthat", repos = "https://cloud.r-project.org")'
fi

echo "==> Installing seurust..."
sed -i 's/\r$//' seurust/configure seurust/cleanup seurust/DESCRIPTION seurust/src/entrypoint.c 2>/dev/null || true
cd seurust
Rscript tools/config.R
cd /workspace
R CMD INSTALL --preclean seurust

echo "==> Rust unit tests..."
cargo test --manifest-path seurust/src/rust/Cargo.toml --lib

echo "==> seurust testthat..."
Rscript -e 'testthat::test_dir("seurust/tests/testthat")'

echo "==> Done."
