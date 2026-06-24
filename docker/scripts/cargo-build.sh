#!/usr/bin/env bash
set -euo pipefail
export R_HOME="$(R RHOME)"
cd /workspace/seurust/src/rust
cargo build --release --lib "$@"
