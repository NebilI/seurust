# Seurat Standard Analysis Benchmark

Times the executable scripts from [ajtimon/seurat-standard-analysis](https://github.com/ajtimon/seurat-standard-analysis) unchanged, comparing:

- **Seurat C++** — default Seurat Rcpp kernels
- **RSeurat Rust** — same scripts after patching Seurat's native entry points to RSeurat

Scripts timed:

| Script | Description |
|--------|-------------|
| `01_pbmc_satija_tutorial.R` | PBMC 3K Satija tutorial |
| `02_gbm_seurat_adapted.R` | GBM object + Harmony integration |
| `03_gbmap_exploration.R` | GBMap atlas exploration |

The runner clones the upstream repo into `upstream/` on first run (or uses an existing clone). Script outputs land in the upstream tree's `output/` and `plots/` directories and are overwritten between backend passes.

## Run With Dev Container

From the repository root:

```sh
docker compose -f docker/docker-compose.yml build rust-dev
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  Rscript benchmarks/seurat_standard_analysis_github/run_benchmark.R
```

## Run As Benchmark Image

```sh
docker compose -f docker/docker-compose.yml build rust-dev
docker build \
  -f benchmarks/seurat_standard_analysis_github/Dockerfile \
  -t rust-seurat-standard-analysis-benchmark .
docker run --rm rust-seurat-standard-analysis-benchmark
```

Mount the working tree for faster iteration:

```sh
docker run --rm -v "${PWD}:/workspace" -w /workspace \
  rust-seurat-standard-analysis-benchmark
```

## Outputs

Timing summary is printed to stdout. Full results are saved to `output/script_timing_results.rds`.

The timing table reports `Seurat / RSeurat` speedup; values above `1.0x` mean the Rust-backed run was faster.
