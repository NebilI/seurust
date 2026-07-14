#!/usr/bin/env bash
# Regenerate vendored Rust crates for CRAN offline builds.
# Run from the repository root: ./seurust/tools/vendor.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUST_DIR="$ROOT/src/rust"

cd "$RUST_DIR"
rm -rf vendor
cargo vendor --locked > vendor-config.toml.tmp
# cargo vendor prints the config on stdout; keep a stable file.
cat > vendor-config.toml <<'EOF'
[source.crates-io]
replace-with = "vendored-sources"

[source.vendored-sources]
directory = "vendor"
EOF
rm -f vendor-config.toml.tmp
rm -f vendor.tar.xz
tar -cJ --no-xattrs -f vendor.tar.xz vendor
rm -rf vendor
ls -lh vendor.tar.xz vendor-config.toml
echo "Vendoring complete. Commit vendor.tar.xz and vendor-config.toml before a CRAN release."
