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
