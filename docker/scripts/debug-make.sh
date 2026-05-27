#!/usr/bin/env bash
set -euo pipefail
cd /workspace
export NOT_CRAN=1
export SEURAT_KEEP_RUST_TARGET=1
Rscript tools/config.R
cd src
R --vanilla -e '
rhome <- R.home()
cat("R_HOME=", rhome, "\n", sep="")
arch <- .Platform$r_arch
makeconf <- file.path(rhome, "etc", arch, "Makeconf")
cat("Makeconf exists:", file.exists(makeconf), makeconf, "\n")
writeLines(sprintf("include %s\ninclude Makevars\nprint-vars:\n\t@echo SHLIB=$(SHLIB)\n\t@echo OBJECTS=$(OBJECTS)\n\t@echo STATLIB=$(STATLIB)", shQuote(makeconf)), "/tmp/test.mk")
system("make -f /tmp/test.mk print-vars")
'
