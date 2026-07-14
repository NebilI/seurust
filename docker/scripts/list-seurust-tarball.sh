#!/usr/bin/env bash
# List seurust source tarball contents (must already exist at repo root).
set -euo pipefail
cd /workspace
TARBALL=$(ls -1t seurust_*.tar.gz 2>/dev/null | head -n 1 || true)
if [ -z "${TARBALL}" ]; then
  echo "No seurust_*.tar.gz found. Run seurust-cran first."
  exit 1
fi
echo "Tarball: ${TARBALL}"
ls -lh "${TARBALL}"
echo ""
echo "=== Top-level entries ==="
tar -tzf "${TARBALL}" | awk -F/ 'NF<=2 {print}' | sort -u
echo ""
echo "=== Full file count ==="
tar -tzf "${TARBALL}" | wc -l
echo ""
echo "=== Unexpected paths (should be empty) ==="
tar -tzf "${TARBALL}" | grep -E 'benchmarks/|docker/|examples/|vignettes/|\.github/|CRAN\.md|vendor\.sh|RcppExports|tests/testthat/test_rust' || echo "(none)"
