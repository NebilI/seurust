# Releasing seurust

One standardized path publishes **GitHub Releases**, **crates.io**, and
**r-universe**. CRAN is not part of this flow.

## What a release is

The version in [`seurust/DESCRIPTION`](seurust/DESCRIPTION) is the source of
truth. The matching Rust crate version lives in
[`seurust/src/rust/Cargo.toml`](seurust/src/rust/Cargo.toml) (and
`Cargo.lock`). The human-readable changelog is
[`seurust/NEWS.md`](seurust/NEWS.md).

When that DESCRIPTION version is new on `main`,
[`.github/workflows/release-seurust.yaml`](.github/workflows/release-seurust.yaml):

1. Creates GitHub Release tag `v<Version>` with notes from the matching
   `NEWS.md` section.
2. Publishes the Rust crate to crates.io.
3. Builds `seurust_<Version>.tar.gz`, attaches it to the release, and upserts
   the seurust entry on `NebilI/NebilI.r-universe.dev`.

Secrets and registry details: [`r-universe/README.md`](r-universe/README.md).

## Prepare a release

```sh
# 1. Draft the top "# seurust X.Y.Z" section in seurust/NEWS.md
#    (improvements, fixes, notes). Keep the heading format exact.

# 2. Bump DESCRIPTION + Cargo.toml + Cargo.lock and verify notes extract:
scripts/prepare-release.sh X.Y.Z

# 3. Optional dry-run of notes only:
scripts/prepare-release.sh X.Y.Z --check-only
# or:
scripts/extract-release-notes.sh X.Y.Z
```

`prepare-release.sh` fails if `NEWS.md` has no section for that version, so
release notes cannot be skipped by accident.

## Ship it

1. Open a PR with the version bump and NEWS.
2. Merge to `main`. A push that changes `seurust/DESCRIPTION` runs
   **Release seurust**.
3. To republish the same version (crate already uploaded is skipped; release
   tag is reused), use **Actions → Release seurust → Run workflow**.

## Versioning

- **Patch** (`0.2.1`): bug fixes only, no kernel behavior changes users rely on.
- **Minor** (`0.3.0`): performance work, new kernels, or packaging changes that
  stay API-compatible with the current Seurat target.
- **Major** (`1.0.0`): breaking R/Rust API or Seurat target changes.

Always keep DESCRIPTION and the Rust crate version identical.

## Helpers

| Script | Role |
|--------|------|
| [`scripts/prepare-release.sh`](scripts/prepare-release.sh) | Bump versions + require NEWS |
| [`scripts/extract-release-notes.sh`](scripts/extract-release-notes.sh) | Print NEWS section for GitHub |
| [`scripts/sync-r-universe.sh`](scripts/sync-r-universe.sh) | Upsert seurust only; preserve harmorust and neighbors |
| [`scripts/test-sync-r-universe.sh`](scripts/test-sync-r-universe.sh) | Regression: harmorust must survive a seurust bump |
| [`scripts/prune-seurust-releases.sh`](scripts/prune-seurust-releases.sh) | Clean stale `v0.1.*` tags if needed |
