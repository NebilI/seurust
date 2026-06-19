# Seurat Standard Analysis Benchmark

This benchmark runs the PBMC 3K Satija tutorial from [ajtimon/seurat-standard-analysis](https://github.com/ajtimon/seurat-standard-analysis) (`code/01_pbmc_satija_tutorial.R`) and compares:

- `Seurat:::` native C++/Rcpp kernels
- `RSeurat::` native Rust/extendr kernels

Upstream script: <https://github.com/ajtimon/seurat-standard-analysis/blob/master/code/01_pbmc_satija_tutorial.R>

Data source: <https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>

The benchmark follows the upstream workflow shape: QC/filter cells, log-normalize, VST variable features, scale all genes, PCA, SNN graph, clustering at resolution 0.5, and UMAP. High-level Seurat object steps stay the same; only ported native kernel calls are swapped between backends.

## Run With Existing Dev Container

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

For faster iteration against your working tree, mount the repo:

```sh
docker run --rm -v "${PWD}:/workspace" -w /workspace \
  rust-seurat-standard-analysis-benchmark
```

On Windows PowerShell, `${PWD}` works for the bind mount.

## Outputs

The script downloads PBMC 3K once into `data/` and writes results into `output/`. Both directories are ignored by git.

Parity checks include cell/gene counts, variable-feature digest, normalized data digest, rounded scaled-data digest, SNN digest, cluster digest, cluster sizes, and SNN non-zero count. The exact scaled-data digest and UMAP digest are reported as informational.

The timing table reports `Seurat / RSeurat`, so values above `1.0x` mean the Rust-backed step was faster.
