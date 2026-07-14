#!/usr/bin/env bash
# Build the CRAN tarball (via check script) and optionally submit to CRAN.
#
# Dry-run (default): build + check only
#   docker compose -f docker/docker-compose.yml run --rm seurust-cran-submit
#
# Submit (uploads tarball; CRAN still emails the maintainer for confirmation):
#   docker compose -f docker/docker-compose.yml run --rm \
#     -e SUBMIT_CRAN=yes \
#     seurust-cran-submit
set -euo pipefail

cd /workspace
sed -i 's/\r$//' docker/scripts/check-seurust-cran.sh 2>/dev/null || true

bash docker/scripts/check-seurust-cran.sh

TARBALL=$(ls -1 seurust_*.tar.gz | head -n 1)
if [ -z "${TARBALL}" ]; then
  echo "No seurust_*.tar.gz found after check."
  exit 1
fi

echo ""
echo "================================================================"
echo " CRAN package (after acceptance):"
echo "   https://cran.r-project.org/package=seurust"
echo " Install once accepted:"
echo "   install.packages(\"seurust\")"
echo " Until then (r-universe):"
echo "   install.packages(\"seurust\","
echo "     repos = c(\"https://NebilI.r-universe.dev\", \"https://cloud.r-project.org\"))"
echo "================================================================"

if [ "${SUBMIT_CRAN:-}" != "yes" ]; then
  echo ""
  echo "Dry-run only. Tarball ready: /workspace/${TARBALL}"
  echo "To upload to CRAN from Docker:"
  echo "  docker compose -f docker/docker-compose.yml run --rm -e SUBMIT_CRAN=yes seurust-cran-submit"
  echo "Or upload manually at https://cran.r-project.org/submit.html"
  echo "  - package: ${TARBALL}"
  echo "  - comments: seurust/cran-comments.md"
  exit 0
fi

echo "==> Installing curl + devtools for CRAN upload..."
Rscript -e 'pkgs <- c("curl", "devtools");
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)];
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")'

echo "==> Attempting CRAN upload of ${TARBALL}..."
# Prefer direct incoming upload used by devtools when prompts cannot be answered
# in non-interactive Docker sessions.
Rscript - <<'EOF'
tarball <- Sys.glob("seurust_*.tar.gz")[[1]]
comments <- paste(readLines("seurust/cran-comments.md"), collapse = "\n")
maintainer <- "Nebil Ibrahim <nebilibrahim@microsoft.com>"

# Mirror the upload endpoint used by remotes/devtools submit helpers.
# If this fails (endpoint/policy change), fall back to the web form instructions.
upload_url <- "https://xmpalantir.wu.ac.at/cransubmit/index2.php"
message("Uploading ", tarball, " to ", upload_url)
h <- curl::handle_setopt(curl::new_handle(), post = TRUE)
curl::handle_setform(
  h,
  uploaded_file = curl::form_file(tarball, type = "application/gzip"),
  name = "Nebil Ibrahim",
  email = "nebilibrahim@microsoft.com",
  comments = comments
)
res <- curl::curl_fetch_memory(upload_url, handle = h)
message("HTTP status: ", res$status_code)
cat(rawToChar(res$content), "\n")
if (res$status_code >= 400) {
  stop("CRAN upload failed with HTTP ", res$status_code,
       ". Upload manually at https://cran.r-project.org/submit.html")
}
message("Upload request completed. Watch email for CRAN confirmation.")
EOF

echo ""
echo "==> Next steps (required):"
echo "  1. Reply to the CRAN confirmation email sent to the Maintainer address."
echo "  2. Address any reviewer follow-ups."
echo "  3. When accepted, users install with: install.packages(\"seurust\")"
echo "  4. Package page: https://cran.r-project.org/package=seurust"
