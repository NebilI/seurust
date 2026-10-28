# r-universe distribution

Publish **seurust** on [r-universe](https://r-universe.dev) so users can install with
`install.packages("seurust", repos = "https://NebilI.r-universe.dev")`.

## One-time setup

1. Create a public GitHub repository named **`NebilI.r-universe.dev`**.
2. Add a root `packages.json` if the registry does not have one yet (see `packages.json.example` in this folder). That example is only the seurust object. Do not copy it over a registry that already lists other packages.
3. Push. r-universe picks up the registry within a few minutes.
4. Dashboard: https://NebilI.r-universe.dev/seurust

## Install for users

```r
install.packages(
  "seurust",
  repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
)
```

## Automated publishing

When you publish a [GitHub Release](https://github.com/NebilI/seurust/releases),
or run **Actions → Release seurust → Run workflow**,
[`publish-seurust-r.yaml`](../.github/workflows/publish-seurust-r.yaml):

1. Builds an R source tarball and attaches it to the release.
2. Upserts the `seurust` object in `packages.json` on `NebilI/NebilI.r-universe.dev` so `branch` is the release tag. Other entries stay, in order.

`NebilI.r-universe.dev` is shared by every package in this universe. [`scripts/sync-r-universe.sh`](../scripts/sync-r-universe.sh) fetches the current array, replaces only the object whose `package` is `seurust` (`url` `https://github.com/NebilI/seurust`, `subdir` `seurust`), and keeps everything else. It creates a one-element array only when `packages.json` is missing (HTTP 404). Any other fetch failure, or a file that is not a JSON array, fails the job and does not commit. A run that does not change the seurust entry does not commit.

CRAN-oriented offline tarballs are built by
[`build-seurust-cran.yaml`](../.github/workflows/build-seurust-cran.yaml).

### One-time GitHub secrets (in `NebilI/seurust`)

| Secret | Used by | Purpose |
|--------|---------|---------|
| `CRATES_IO_TOKEN` | `publish-seurust-crate.yaml` | Publish the `seurust` crate to [crates.io](https://crates.io) |
| `R_UNIVERSE_REGISTRY_TOKEN` | `publish-seurust-r.yaml` | PAT with `repo` scope to upsert the seurust entry in `NebilI/NebilI.r-universe.dev` `packages.json` |

Create the crates.io token at https://crates.io/settings/tokens (needs `publish-new` /
`publish-update` for the `seurust` crate).

### Release checklist

1. If Rust deps changed, refresh vendoring via Docker:
   `docker compose -f docker/docker-compose.yml run --rm rust-dev bash docker/scripts/vendor-seurust.sh`
2. Bump `Version` in [`seurust/DESCRIPTION`](../seurust/DESCRIPTION) and `version` in
   [`seurust/src/rust/Cargo.toml`](../seurust/src/rust/Cargo.toml) together.
3. Update [`seurust/NEWS.md`](../seurust/NEWS.md) and [`seurust/cran-comments.md`](../seurust/cran-comments.md).
4. Validate CRAN tarball:
   `docker compose -f docker/docker-compose.yml run --rm seurust-cran`
5. Merge a version bump to `main`. That tags the release. To republish the current version without a new commit, use **Actions → Release seurust → Run workflow**.
6. That run updates r-universe and crates.io. Trigger `build-seurust-cran` when you also want a CI CRAN tarball artifact.
7. For CRAN: use the Docker-built `seurust_*.tar.gz`, then `devtools::submit_cran()` from `seurust/`
   (or upload via the CRAN web form).
