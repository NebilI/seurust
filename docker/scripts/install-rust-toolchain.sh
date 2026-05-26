#!/usr/bin/env bash
set -euo pipefail

export RUSTUP_HOME="${RUSTUP_HOME:-/usr/local/rustup}"
export CARGO_HOME="${CARGO_HOME:-/usr/local/cargo}"
export PATH="${CARGO_HOME}/bin:${PATH}"

if ! command -v rustup >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --no-modify-path --profile minimal
fi

rustup default stable
rustup component add rustfmt clippy
rustup target add x86_64-unknown-linux-gnu 2>/dev/null || true

rustc --version
cargo --version
