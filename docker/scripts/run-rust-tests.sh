#!/usr/bin/env bash
set -euo pipefail
cd /workspace
sed -i 's/\r$//' docker/scripts/build-and-test-rust.sh
export NOT_CRAN=1
export SEURAT_KEEP_RUST_TARGET=1
Rscript docker/scripts/install-imports.R
bash docker/scripts/build-and-test-rust.sh
