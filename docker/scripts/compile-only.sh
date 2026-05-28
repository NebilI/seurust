#!/usr/bin/env bash
set -euo pipefail
cd /workspace
export NOT_CRAN=1
Rscript docker/scripts/install-imports.R
Rscript -e 'pkgbuild::compile_dll(debug = FALSE, compile_attributes = FALSE)'
Rscript -e 'devtools::load_all(recompile = FALSE)'
echo "BUILD OK (Seurat C++)"
