#!/usr/bin/env bash
# Build the CRAN tarball (via check script) and optionally submit to CRAN.
#
# Dry-run (default): build + check only
#   docker compose -f docker/docker-compose.yml run --rm seurust-cran-submit
#
# Submit (two-step CRAN web upload; maintainer must confirm email):
#   docker compose -f docker/docker-compose.yml run --rm \
#     -e SUBMIT_CRAN=yes \
#     seurust-cran-submit
#
# Upload an already-built tarball without re-checking:
#   docker compose -f docker/docker-compose.yml run --rm \
#     -e SUBMIT_CRAN=yes -e SKIP_CHECK=yes \
#     seurust-cran-submit
set -euo pipefail

cd /workspace
sed -i 's/\r$//' docker/scripts/check-seurust-cran.sh 2>/dev/null || true

if [ "${SKIP_CHECK:-}" = "yes" ]; then
  echo "==> SKIP_CHECK=yes; using existing seurust_*.tar.gz"
else
  bash docker/scripts/check-seurust-cran.sh
fi

TARBALL=$(ls -1t seurust_*.tar.gz 2>/dev/null | head -n 1 || true)
if [ -z "${TARBALL}" ]; then
  echo "No seurust_*.tar.gz found."
  exit 1
fi

echo ""
echo "================================================================"
echo " CRAN package page (after acceptance):"
echo "   https://cran.r-project.org/package=seurust"
echo " Install once accepted:"
echo "   install.packages(\"seurust\")"
echo " Until then (r-universe):"
echo "   install.packages(\"seurust\","
echo "     repos = c(\"https://NebilI.r-universe.dev\", \"https://cloud.r-project.org\"))"
echo " Tarball: /workspace/${TARBALL}"
echo "================================================================"

if [ "${SUBMIT_CRAN:-}" != "yes" ]; then
  echo ""
  echo "Dry-run only."
  echo "To upload to CRAN from Docker:"
  echo "  docker compose -f docker/docker-compose.yml run --rm -e SUBMIT_CRAN=yes seurust-cran-submit"
  echo "Or upload manually at https://cran.r-project.org/submit.html"
  echo "  - package: ${TARBALL}"
  echo "  - comments: seurust/cran-comments.md"
  exit 0
fi

echo "==> Installing curl + httr2 for CRAN upload..."
Rscript -e 'pkgs <- c("curl", "httr2");
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)];
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")'

echo "==> Uploading ${TARBALL} to CRAN (devtools two-step flow)..."
export SEURUST_TARBALL="${TARBALL}"
Rscript - <<'EOF'
tarball <- Sys.getenv("SEURUST_TARBALL")
stopifnot(nzchar(tarball), file.exists(tarball))
comments <- paste(readLines("seurust/cran-comments.md", warn = FALSE), collapse = "\n")
desc <- read.dcf("seurust/DESCRIPTION")
# Prefer Authors@R cre email when present.
authors <- desc[1, "Authors@R"]
email <- sub('.*email\\s*=\\s*"([^"]+)".*', "\\1", authors)
if (!grepl("@", email)) {
  email <- "nebil080298@gmail.com"
}
name <- "Nebil Ibrahim"
upload_url <- "https://xmpalantir.wu.ac.at/cransubmit/index2.php"

message("Maintainer email for CRAN: ", email)
message("Step 1/2: upload ", tarball)
req <- httr2::request(upload_url)
req <- httr2::req_body_multipart(
  req,
  pkg_id = "",
  name = name,
  email = email,
  uploaded_file = curl::form_file(tarball, type = "application/x-gzip"),
  comment = comments,
  upload = "Upload package"
)
resp <- httr2::req_perform(req)
final_url <- httr2::resp_url(resp)
message("Redirect/URL: ", final_url)
parsed <- httr2::url_parse(final_url)
pkg_id <- parsed$query$pkg_id
if (is.null(pkg_id) || !nzchar(pkg_id)) {
  cat(httr2::resp_body_string(resp), "\n")
  stop("CRAN upload did not return pkg_id. Upload manually at https://cran.r-project.org/submit.html")
}

message("Step 2/2: confirm submission (pkg_id=", pkg_id, ")")
req2 <- httr2::request(upload_url)
req2 <- httr2::req_body_multipart(
  req2,
  pkg_id = pkg_id,
  name = name,
  email = email,
  policy_check = "1/",
  submit = "Submit package"
)
resp2 <- httr2::req_perform(req2)
final2 <- httr2::resp_url(resp2)
message("Confirm URL: ", final2)
parsed2 <- httr2::url_parse(final2)
if (is.null(parsed2$query$submit) || parsed2$query$submit != "1") {
  cat(httr2::resp_body_string(resp2), "\n")
  stop("CRAN confirmation step failed. Try the web form: https://cran.r-project.org/submit.html")
}

writeLines(
  c(
    paste0("version: ", sub("^seurust_(.*)\\.tar\\.gz$", "\\1", basename(tarball))),
    paste0("tarball: ", tarball),
    paste0("submitted_at: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("pkg_id: ", pkg_id),
    paste0("email: ", email)
  ),
  "seurust/CRAN-SUBMISSION"
)
message("Package submission successful.")
message("Check ", email, " for the CRAN confirmation link, then reply/confirm.")
message("After acceptance: https://cran.r-project.org/package=seurust")
EOF

echo ""
echo "==> Next steps (required):"
echo "  1. Open the confirmation link in the CRAN email (Gmail: nebil080298@gmail.com)."
echo "  2. Do NOT confirm any older CRAN email sent to a Microsoft address."
echo "  3. Address any reviewer follow-ups."
echo "  4. When accepted: install.packages(\"seurust\")"
echo "  5. Package page: https://cran.r-project.org/package=seurust"
