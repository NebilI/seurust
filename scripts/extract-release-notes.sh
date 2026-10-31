#!/usr/bin/env bash
# Extract the NEWS.md section for one seurust version.
#
# Usage:
#   scripts/extract-release-notes.sh [VERSION]
#   scripts/extract-release-notes.sh 0.3.0
#
# When VERSION is omitted, reads Version from seurust/DESCRIPTION.
# Prints Markdown suitable for a GitHub Release body. Exits non-zero if the
# section is missing or empty.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
news="$root/seurust/NEWS.md"
desc="$root/seurust/DESCRIPTION"

version="${1:-}"
if [[ -z "$version" ]]; then
  version="$(awk '/^Version:[[:space:]]*/ { print $2; exit }' "$desc")"
fi
if [[ -z "$version" ]]; then
  echo "Could not determine version (pass it or set seurust/DESCRIPTION)." >&2
  exit 1
fi
if [[ ! -f "$news" ]]; then
  echo "Missing $news" >&2
  exit 1
fi

# Match "# seurust X.Y.Z" (optional trailing whitespace) through the next
# same-level heading or EOF. Keep the heading so the release body is self-
# contained on GitHub.
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT
awk -v ver="$version" '
  BEGIN { want = "# seurust " ver }
  $0 ~ /^# seurust / {
    if (capturing) { exit }
    if ($0 == want || $0 ~ ("^# seurust " ver "[[:space:]]*$")) {
      capturing = 1
    }
  }
  capturing { print }
  END {
    if (!capturing) {
      print "No NEWS.md section for seurust " ver > "/dev/stderr"
      exit 1
    }
  }
' "$news" >"$tmp"
if ! grep -q '[^[:space:]]' "$tmp"; then
  echo "NEWS.md section for seurust ${version} is empty" >&2
  exit 1
fi
cat "$tmp"
