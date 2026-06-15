# PBMC 3K Seurat GitHub Benchmark

This benchmark runs a Docker-friendly version of the public Satija Lab PBMC 3K Seurat tutorial and compares:

- `Seurat:::` native C++/Rcpp kernels
- `RSeurat::` native Rust/extendr kernels

Upstream workflow: <https://github.com/satijalab/seurat/blob/HEAD/vignettes/pbmc3k_tutorial.Rmd>

Data source: <https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>

The script follows the tutorial shape: create the object, QC/filter cells, log-normalize, select variable features, scale, run PCA, build an SNN graph, cluster, and run UMAP. Since `RSeurat` is currently a companion package for ported native kernels, the high-level Seurat object workflow stays the same and only the backend calls for ported kernels are swapped.

## Run With Existing Dev Container

From the repository root:

```sh
docker compose -f docker/docker-compose.yml build rust-dev
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  Rscript benchmarks/pbmc3k_seurat_github/run_benchmark.R
```

## Run As Benchmark Image

Build the existing Rust dev base image first, then build this benchmark image:

```sh
docker compose -f docker/docker-compose.yml build rust-dev
docker build \
  -f benchmarks/pbmc3k_seurat_github/Dockerfile \
  -t rust-seurat-pbmc3k-benchmark .
docker run --rm rust-seurat-pbmc3k-benchmark
```

For faster iteration against your working tree, mount the repo instead of using the copied checkout baked into the image:

```sh
docker run --rm -v "${PWD}:/workspace" -w /workspace \
  rust-seurat-pbmc3k-benchmark
```

On Windows PowerShell, `${PWD}` works for the bind mount.

## Outputs

The script downloads PBMC 3K once into `data/` and writes results into `output/`. Both directories are intentionally ignored by git.

Parity checks include cell/gene counts, variable-feature digest, normalized data digest, rounded scaled-data digest, SNN digest, cluster digest, cluster sizes, and SNN non-zero count. The exact scaled-data digest is still reported as informational because tiny native floating-point differences can appear there while preserving the downstream graph and clustering outputs. UMAP digest is also informational because it can vary across `uwot`/threading/platform details even with a fixed seed.

The timing table reports `Seurat / RSeurat`, so values above `1.0x` mean the Rust-backed step was faster.
