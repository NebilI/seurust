#!/usr/bin/env bash
# Refresh vendored Rust crates for seurust CRAN builds inside rust-dev.
#
# From repo root:
#   docker compose -f docker/docker-compose.yml run --rm \
#     rust-dev bash docker/scripts/vendor-seurust.sh
set -euo pipefail

cd /workspace
sed -i 's/\r$//' seurust/tools/vendor.sh 2>/dev/null || true
chmod +x seurust/tools/vendor.sh
bash seurust/tools/vendor.sh
