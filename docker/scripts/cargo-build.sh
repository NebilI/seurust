#!/usr/bin/env bash
set -euo pipefail
export R_HOME="$(R RHOME)"
cd /workspace/RSeurat/src/rust
cargo build --release --lib "$@"
