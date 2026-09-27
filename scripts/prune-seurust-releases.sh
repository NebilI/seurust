#!/usr/bin/env bash
# Remove stale seurust GitHub release tags (keep the version in seurust/DESCRIPTION).
set -euo pipefail

repo="${GITHUB_REPOSITORY:-NebilI/seurust}"
root="$(cd "$(dirname "$0")/.." && pwd)"
keep="$(awk '/^Version:[[:space:]]*/ { print "v" $2; exit }' "$root/seurust/DESCRIPTION")"

if [[ -z "$keep" ]]; then
  echo "Could not read Version from seurust/DESCRIPTION" >&2
  exit 1
fi

echo "Keeping release tag: $keep"
mapfile -t tags < <(gh api "repos/${repo}/tags?per_page=100" --jq '.[].name' | rg '^v0\.1\.' || true)

for tag in "${tags[@]:-}"; do
  [[ "$tag" == "$keep" ]] && continue
  echo "Deleting release and tag: $tag"
  gh release delete "$tag" -y --repo "$repo" 2>/dev/null || true
  git push origin ":refs/tags/${tag}" 2>/dev/null || true
done

echo "Done. Remaining seurust tags:"
gh api "repos/${repo}/tags?per_page=100" --jq '.[].name' | rg '^v0\.1\.' || echo "(none)"
