#!/usr/bin/env bash
# Standardize a seurust release bump: sync versions, require NEWS, dry-run notes.
#
# Usage:
#   scripts/prepare-release.sh 0.3.0
#   scripts/prepare-release.sh 0.3.0 --check-only
#
# Without --check-only this writes:
#   - seurust/DESCRIPTION Version
#   - seurust/src/rust/Cargo.toml version
#   - seurust/src/rust/Cargo.lock package version for seurust
#
# It always verifies that seurust/NEWS.md has a "# seurust <version>" section
# and prints the GitHub release notes that release-seurust.yaml will use.
#
# After this script (and a NEWS edit if needed):
#   1. Commit the bump on a PR.
#   2. Merge to main (or push DESCRIPTION) — Release seurust publishes.
#   3. Or run Actions → Release seurust to republish the current version.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
desc="$root/seurust/DESCRIPTION"
cargo_toml="$root/seurust/src/rust/Cargo.toml"
cargo_lock="$root/seurust/src/rust/Cargo.lock"
news="$root/seurust/NEWS.md"
extract="$root/scripts/extract-release-notes.sh"

usage() {
  echo "Usage: $0 <version> [--check-only]" >&2
  echo "  version      Semver like 0.3.0 (no leading v)" >&2
  echo "  --check-only Validate NEWS and print notes; do not rewrite files" >&2
}

version=""
check_only=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check-only) check_only=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      if [[ -n "$version" ]]; then
        echo "Unexpected argument: $1" >&2
        usage
        exit 1
      fi
      version="$1"
      shift
      ;;
  esac
done

if [[ -z "$version" ]]; then
  usage
  exit 1
fi
if [[ ! "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-][A-Za-z0-9.+-]+)?$ ]]; then
  echo "Version must look like 0.3.0 (got: $version)" >&2
  exit 1
fi
if [[ "$version" == v* ]]; then
  echo "Pass the bare version without a leading v (got: $version)" >&2
  exit 1
fi

for path in "$desc" "$cargo_toml" "$cargo_lock" "$news" "$extract"; do
  if [[ ! -e "$path" ]]; then
    echo "Missing required path: $path" >&2
    exit 1
  fi
done

current="$(awk '/^Version:[[:space:]]*/ { print $2; exit }' "$desc")"
echo "Current DESCRIPTION version: ${current:-unknown}"
echo "Target release version:      $version"

if [[ "$check_only" -eq 0 ]]; then
  # DESCRIPTION
  awk -v ver="$version" '
    /^Version:[[:space:]]*/ { print "Version: " ver; next }
    { print }
  ' "$desc" >"$desc.tmp"
  mv "$desc.tmp" "$desc"

  # Cargo.toml package version only (first occurrence under [package])
  awk -v ver="$version" '
    BEGIN { in_pkg = 0; done = 0 }
    /^\[package\]/ { in_pkg = 1 }
    /^\[/ && !/^\[package\]/ { in_pkg = 0 }
    in_pkg && !done && /^version[[:space:]]*=/ {
      print "version = \"" ver "\""
      done = 1
      next
    }
    { print }
  ' "$cargo_toml" >"$cargo_toml.tmp"
  mv "$cargo_toml.tmp" "$cargo_toml"

  # Cargo.lock local package
  awk -v ver="$version" '
    BEGIN { hit = 0 }
    /^name = "seurust"$/ { hit = 1; print; next }
    hit && /^version = / {
      print "version = \"" ver "\""
      hit = 0
      next
    }
    { print }
  ' "$cargo_lock" >"$cargo_lock.tmp"
  mv "$cargo_lock.tmp" "$cargo_lock"

  echo "Updated DESCRIPTION, Cargo.toml, and Cargo.lock to $version"
else
  echo "Check-only: leaving version files unchanged"
fi

echo
echo "---- GitHub release notes (from seurust/NEWS.md) ----"
"$extract" "$version"
echo "---- end ----"
echo
echo "Next:"
echo "  1. Confirm NEWS.md lists every improvement and fix for $version."
echo "  2. Commit and open a PR; merge to main to publish."
echo "  3. Or: Actions → Release seurust (reuses tag v${version} if present)."
