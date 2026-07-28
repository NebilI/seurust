# Publishing seurust (r-universe, GitHub Releases, CRAN)

This guide is the canonical process for making **seurust** publicly installable.
All local build/check/submit steps run **inside Docker Compose** — do not use a
host R or Cargo toolchain.

## Quick commands (from repo root)

```sh
# Build images (once)
docker compose -f docker/docker-compose.yml build rcpp-dev rust-dev

# CRAN-style offline build + R CMD check --as-cran
docker compose -f docker/docker-compose.yml run --rm seurust-cran

# Refresh vendored Rust crates (needed when Cargo.lock / deps change)
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  bash docker/scripts/vendor-seurust.sh

# Prepare a CRAN submission (builds tarball; upload only with SUBMIT_CRAN=yes READY_TO_PUBLISH=yes)
docker compose -f docker/docker-compose.yml run --rm \
  -e SUBMIT_CRAN=yes -e READY_TO_PUBLISH=yes \
  seurust-cran-submit
```

## Public install channels

| Channel | When available | User install |
|---------|----------------|--------------|
| **r-universe** | After registry repo is live + first build | `install.packages("seurust", repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org"))` |
| **GitHub Release** | On every GitHub Release | `install.packages("seurust_0.1.0.tar.gz", repos = NULL, type = "source")` |
| **CRAN** | Only after CRAN **accepts** the submission | `install.packages("seurust")` |
| **crates.io** | On release (Rust crate) | mostly for packaging; users install the R package |

### CRAN package page

https://cran.r-project.org/package=seurust

**Status: submitting 0.1.0.** Upload with
`SUBMIT_CRAN=yes READY_TO_PUBLISH=yes`. The CRAN page 404s until acceptance;
use r-universe/GitHub until then.

---

## How CRAN publishing actually works

CRAN is a curated repository. There is **no fully automatic “push to CRAN” API**
like crates.io or npm. Maintainers submit a source tarball; humans review it.

### 1. Prepare a release candidate

1. Bump versions together:
   - `seurust/DESCRIPTION` → `Version: x.y.z` (no `.9000` suffix)
   - `seurust/src/rust/Cargo.toml` → matching `version`
2. Update `seurust/NEWS.md` and `seurust/cran-comments.md`.
3. If Rust dependencies changed, vendor again (Docker command above).
4. Run the CRAN check service until you have **no ERRORs** (WARNINGs/NOTEs must be explained in `cran-comments.md`).

### 2. Submit the tarball

Submission uploads `seurust_x.y.z.tar.gz` to CRAN’s incoming area (via
`devtools::submit_cran()` or the [CRAN web form](https://cran.r-project.org/submit.html)).

Our Docker submit service wraps the check + `devtools::submit_cran()` path:

```sh
docker compose -f docker/docker-compose.yml run --rm \
  -e SUBMIT_CRAN=yes -e READY_TO_PUBLISH=yes \
  seurust-cran-submit
```

### 3. Confirm by email (required)

CRAN emails the **Maintainer** address from `DESCRIPTION`
(`nbi@alumni.princeton.edu`). You must reply to confirm the submission.
Without that reply, the package never enters review.

If an earlier submission used a different maintainer address, **do not confirm that
email**. Confirm only the message sent to `nbi@alumni.princeton.edu`.

### 4. Respond to reviewer feedback

CRAN may ask for changes (docs, portability, size, policies). Fix on a branch,
re-check with Docker, bump a patch version if needed, and resubmit.

### 5. Acceptance and mirrors

After acceptance:

- Package page: https://cran.r-project.org/package=seurust
- Users run: `install.packages("seurust")`
- Mirrors propagate over hours; win/mac binaries appear after CRAN builders run
  (source install always needs Rust + Cargo).

Typical first-submission turnaround is **a few days to a couple of weeks**,
depending on reviewer load and issues found.

---

## CI / release pipelines

| Workflow | Trigger | What it does |
|----------|---------|--------------|
| `seurust_checks.yaml` | Pull requests / pushes touching seurust or Docker packaging | R CMD check + testthat for seurust |
| `build-seurust-cran.yaml` | `workflow_dispatch`, GitHub Release | Docker Compose CRAN build/check; optional CRAN upload |
| `publish-seurust-r.yaml` | GitHub Release | Release tarball + sync `NebilI.r-universe.dev` |
| `publish-seurust-crate.yaml` | GitHub Release | `cargo publish` to crates.io |

### Update / resubmit to CRAN from GitHub Actions

1. Bump `seurust/DESCRIPTION` (and matching `seurust/src/rust/Cargo.toml`) on a PR; merge after `seurust Checks` is green.
2. On `main`: **Actions → Build / submit seurust to CRAN → Run workflow**.
3. Leave **submit_to_cran** unchecked for a dry-run (artifact only), or check it to upload.
4. Confirm the email sent to `nbi@alumni.princeton.edu`, then watch https://cran.r-project.org/package=seurust.

Local equivalent (Docker):

```sh
docker compose -f docker/docker-compose.yml run --rm \
  -e SUBMIT_CRAN=yes -e READY_TO_PUBLISH=yes \
  seurust-cran-submit
```

### One-time secrets (`NebilI/seurust` → Settings → Secrets)

| Secret | Purpose |
|--------|---------|
| `R_UNIVERSE_REGISTRY_TOKEN` | PAT (`repo`) to update `NebilI/NebilI.r-universe.dev` |
| `CRATES_IO_TOKEN` | Publish the Rust crate |

CRAN submission itself uses email confirmation, not a GitHub secret.

### Suggested release flow

1. Land changes on `main` (PR checks via `seurust_checks.yaml`).
2. Tag a release (for example `v0.1.1`) → r-universe + crates.io + CRAN tarball artifact.
3. Run **Build / submit seurust to CRAN** with `submit_to_cran=true` when ready.
4. Confirm the CRAN email and watch https://cran.r-project.org/package=seurust.

---

## Policy notes for this package

- **Rust toolchain**: `SystemRequirements: Cargo …, rustc (>= 1.81)`.
- **Vendoring**: `src/rust/vendor.tar.xz` ships so CRAN builds offline.
- **GNU make**: `.NOTPARALLEL` in Makevars is intentional (extendr template); called out in `cran-comments.md`.
- **Suggests Seurat**: optional parity comparisons; not required to install seurust.
