# r-universe distribution

Use this folder as a template for publishing **RSeurat** on [r-universe](https://r-universe.dev).

## One-time setup

1. Create a public GitHub repository named **`NebilI.r-universe.dev`** (replace `NebilI` with your GitHub username).
2. Copy `packages.json.example` into that repo as **`packages.json`**.
3. Adjust `branch` if you publish from a branch other than `main`.
4. Push the file. r-universe will pick up the registry automatically within a few minutes.

## Install for users

After the universe is live:

```r
install.packages(
  "RSeurat",
  repos = c("https://NebilI.r-universe.dev", "https://cloud.r-project.org")
)
```

Or enable the repository once per session:

```r
options(repos = c(NebilI = "https://NebilI.r-universe.dev", CRAN = "https://cloud.r-project.org"))
install.packages("RSeurat")
```

## Notes

- r-universe builds from source; users still need a Rust toolchain unless you publish pre-built binaries via a custom workflow.
- For experimental branches, set `"branch": "feature/rust-rewrite"` in `packages.json`.
- Dashboard: `https://NebilI.r-universe.dev/RSeurat`

## Automated publishing

When you publish a [GitHub Release](https://github.com/NebilI/Rust-Seurat/releases), the
[`publish-rseurat-r.yaml`](../.github/workflows/publish-rseurat-r.yaml) workflow:

1. Builds an R source tarball and attaches it to the release.
2. Updates `packages.json` in `NebilI.r-universe.dev` to point at the release tag.

### One-time GitHub secrets

Add these under **Settings → Secrets and variables → Actions** in `NebilI/Rust-Seurat`:

| Secret | Used by | Purpose |
|--------|---------|---------|
| `CRATES_IO_TOKEN` | `publish-rseurat-crate.yaml` | Publish the `rseurat` crate to [crates.io](https://crates.io) |
| `R_UNIVERSE_REGISTRY_TOKEN` | `publish-rseurat-r.yaml` | PAT with `repo` scope to push `packages.json` to `NebilI/NebilI.r-universe.dev` |

Create the crates.io token at https://crates.io/settings/tokens (needs `publish-new` /
`publish-update` for the `rseurat` crate).

### Release checklist

1. Bump `Version` in [`RSeurat/DESCRIPTION`](../RSeurat/DESCRIPTION) and `version` in
   [`RSeurat/src/rust/Cargo.toml`](../RSeurat/src/rust/Cargo.toml) together.
2. Tag the release (for example `v0.1.0`) and publish a GitHub Release from that tag.
3. Both publish workflows run automatically; you can also trigger them manually from the
   Actions tab.
