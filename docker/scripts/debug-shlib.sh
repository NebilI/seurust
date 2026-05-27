#!/usr/bin/env bash
set -euo pipefail
cd /workspace
export NOT_CRAN=1
export SEURAT_KEEP_RUST_TARGET=1
Rscript tools/config.R
cd src
echo "=== OBJECTS in dir ==="
ls *.cpp
echo "=== Running R CMD SHLIB ==="
R CMD SHLIB . > /tmp/shlib.log 2>&1 || { cat /tmp/shlib.log; exit 1; }
cat /tmp/shlib.log
echo "=== artifacts ==="
ls -la *.o Seurat.so 2>&1 || true
