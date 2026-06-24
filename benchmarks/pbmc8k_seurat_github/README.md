# PBMC 8K Seurat Benchmark

This benchmark runs a Docker-friendly version of the standard Seurat guided-clustering workflow on the 10x Genomics PBMC 8K dataset and compares:

- `Seurat:::` native C++/Rcpp kernels
- `seurust::` native Rust/extendr kernels

Upstream workflow shape: [Seurat PBMC 3K guided tutorial](https://github.com/satijalab/seurat/blob/HEAD/vignettes/pbmc3k_tutorial.Rmd) (the same analysis steps applied to a larger public dataset).

Data source: [10x Genomics PBMC 8K](https://www.10xgenomics.com/datasets/8-k-pbm-cs-from-a-healthy-donor-2-standard-2-1-0)

Related R projects using the same Seurat + 10x PBMC pattern:

- [Bioconductor TENxPBMCData](https://bioconductor.org/packages/TENxPBMCData) (includes `pbmc8k` among other 10x PBMC sets)
- Community Seurat pipelines such as [NoWon1/scRNA_seq_analysis](https://github.com/NoWon1/scRNA_seq_analysis)

The script follows the tutorial shape: create the object, QC/filter cells, log-normalize, select variable features, scale, run PCA, build an SNN graph, cluster, and run UMAP. Since `seurust` is currently a companion package for ported native kernels, the high-level Seurat object workflow stays the same and only the backend calls for ported kernels are swapped.

## Run With Existing Dev Container

From the repository root:

```sh
docker compose -f docker/docker-compose.yml build rust-dev
docker compose -f docker/docker-compose.yml run --rm rust-dev \
  Rscript benchmarks/pbmc8k_seurat_github/run_benchmark.R
```

## Run As Benchmark Image

Build the existing Rust dev base image first, then build this benchmark image:

```sh
docker compose -f docker/docker-compose.yml build rust-dev
docker build \
  -f benchmarks/pbmc8k_seurat_github/Dockerfile \
  -t seurust-pbmc8k-benchmark .
docker run --rm seurust-pbmc8k-benchmark
```

For faster iteration against your working tree, mount the repo instead of using the copied checkout baked into the image:

```sh
docker run --rm -v "${PWD}:/workspace" -w /workspace \
  seurust-pbmc8k-benchmark
```

On Windows PowerShell, `${PWD}` works for the bind mount.

## Outputs

The script downloads PBMC 8K once into `data/` and writes results into `output/`. Both directories are intentionally ignored by git.

Parity checks include cell/gene counts, variable-feature digest, normalized data digest, rounded scaled-data digest, SNN digest, cluster digest, cluster sizes, and SNN non-zero count. The exact scaled-data digest is still reported as informational because tiny native floating-point differences can appear there while preserving the downstream graph and clustering outputs. UMAP digest is also informational because it can vary across `uwot`/threading/platform details even with a fixed seed.

The timing table reports `Seurat / seurust`, so values above `1.0x` mean the Rust-backed step was faster.

## Notes

- PBMC 8K uses the GRCh38 reference (PBMC 3K uses hg19), so this benchmark stresses a larger matrix with a different gene annotation.
- For an integration-style workflow (e.g. Seurat's IFN-beta `ifnb` tutorial via `SeuratData`), see the integration introduction vignette; that dataset requires the SeuratData repository, which may not be reachable from all CI/Docker environments.
