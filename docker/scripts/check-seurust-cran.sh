#!/usr/bin/env bash
# Build and R CMD check seurust in CRAN-offline mode inside rust-dev.
#
# From repo root:
#   docker compose -f docker/docker-compose.yml run --rm \
#     -e NOT_CRAN= -e SEURAT_KEEP_RUST_TARGET= \
#     rust-dev bash docker/scripts/check-seurust-cran.sh
set -euo pipefail

cd /workspace

# CRAN mode: unset NOT_CRAN so tools/config.R emits --offline cargo flags.
unset NOT_CRAN || true
unset SEURAT_KEEP_RUST_TARGET || true
export PATH="${PATH}:${HOME}/.cargo/bin:/usr/local/cargo/bin"

echo "==> Ensuring check tooling (checkbashisms)..."
if ! command -v checkbashisms >/dev/null 2>&1; then
  apt-get update -qq
  apt-get install -y --no-install-recommends devscripts >/dev/null
fi

echo "==> Ensuring R packages for check/vignettes..."
Rscript -e 'pkgs <- c("Matrix", "RcppEigen", "knitr", "rmarkdown", "testthat", "Seurat");
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)];
  if (length(miss)) {
    message("Installing: ", paste(miss, collapse = ", "));
    install.packages(miss, repos = c("https://cloud.r-project.org", "https://satijalab.r-universe.dev"));
  }'

echo "==> Normalizing line endings for configure/scripts..."
sed -i 's/\r$//' \
  seurust/configure \
  seurust/cleanup \
  seurust/tools/config.R \
  seurust/tools/msrv.R \
  seurust/tools/vendor.sh \
  2>/dev/null || true
chmod +x seurust/configure seurust/cleanup seurust/tools/vendor.sh

echo "==> Ensuring vendored crates (CRAN offline)..."
if [ ! -f seurust/src/rust/vendor.tar.xz ] || [ ! -f seurust/src/rust/vendor-config.toml ]; then
  bash seurust/tools/vendor.sh
else
  echo "    using existing seurust/src/rust/vendor.tar.xz"
fi

echo "==> Cleaning previous build/check artifacts..."
rm -rf seurust.Rcheck seurust_*.tar.gz

echo "==> Configuring seurust for CRAN..."
cd seurust
Rscript tools/config.R
cd /workspace

echo "==> R CMD build seurust..."
rm -f seurust_*.tar.gz
R CMD build seurust
TARBALL=$(ls -1 seurust_*.tar.gz | head -n 1)
echo "    built ${TARBALL}"

echo "==> R CMD check --as-cran ${TARBALL}..."
# Keep the check directory for inspection on failure.
set +e
R CMD check --as-cran --no-manual "${TARBALL}"
status=$?
set -e

if [ -f seurust.Rcheck/00check.log ]; then
  echo "==> check summary:"
  tail -n 5 seurust.Rcheck/00check.log
fi

if [ "${status}" -ne 0 ]; then
  echo "==> check failed; last log lines:"
  if [ -f seurust.Rcheck/00check.log ]; then
    tail -n 80 seurust.Rcheck/00check.log
  fi
  exit "${status}"
fi

# Treat ERROR as failure even if exit code is odd; also surface WARNINGs.
if grep -q '^Status:.*ERROR' seurust.Rcheck/00check.log; then
  echo "==> ERROR found in check log"
  exit 1
fi

echo "==> CRAN check finished (see WARNINGs/NOTEs above if any)."
ls -lh "${TARBALL}"
echo "Tarball ready for submission: /workspace/${TARBALL}"
