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

# Prepare a CRAN submission (builds tarball; upload only with SUBMIT_CRAN=yes)
docker compose -f docker/docker-compose.yml run --rm \
  -e SUBMIT_CRAN=yes \
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

`seurust` **0.1.0** has been submitted to CRAN (maintainer:
`nebil080298@gmail.com`). The package page and `install.packages("seurust")`
become available after CRAN acceptance and mirror sync.

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
  -e SUBMIT_CRAN=yes \
  seurust-cran-submit
```

### 3. Confirm by email (required)

CRAN emails the **Maintainer** address from `DESCRIPTION`
(`nebil080298@gmail.com`). You must reply to confirm the submission.
Without that reply, the package never enters review.

If an earlier submission used a different maintainer address, **do not confirm that
email**. Resubmit with the correct `Authors@R` email (via Docker
`seurust-cran-submit`) and confirm only the Gmail message.

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
| `build-seurust-cran.yaml` | `workflow_dispatch`, GitHub Release | Docker Compose CRAN build/check; uploads `seurust_*.tar.gz` artifact |
| `publish-seurust-r.yaml` | GitHub Release | Release tarball + sync `NebilI.r-universe.dev` |
| `publish-seurust-crate.yaml` | GitHub Release | `cargo publish` to crates.io |

### One-time secrets (`NebilI/seurust` → Settings → Secrets)

| Secret | Purpose |
|--------|---------|
| `R_UNIVERSE_REGISTRY_TOKEN` | PAT (`repo`) to update `NebilI/NebilI.r-universe.dev` |
| `CRATES_IO_TOKEN` | Publish the Rust crate |

CRAN submission itself uses email confirmation, not a GitHub secret.

### Suggested release flow

1. Land changes on `main`.
2. Tag `v0.1.0` and publish a GitHub Release → r-universe + crates.io + artifacts.
3. Run `seurust-cran` / `seurust-cran-submit` when ready for CRAN.
4. Confirm the CRAN email and watch https://cran.r-project.org/package=seurust.

---

## Policy notes for this package

- **Rust toolchain**: `SystemRequirements: Cargo …, rustc (>= 1.81)`.
- **Vendoring**: `src/rust/vendor.tar.xz` ships so CRAN builds offline.
- **GNU make**: `.NOTPARALLEL` in Makevars is intentional (extendr template); called out in `cran-comments.md`.
- **Suggests Seurat**: optional parity comparisons; not required to install seurust.
