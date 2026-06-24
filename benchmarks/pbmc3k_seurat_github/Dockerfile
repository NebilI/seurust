ARG BASE_IMAGE=seurust-rust:dev
FROM ${BASE_IMAGE}

WORKDIR /workspace
COPY . /workspace

CMD ["Rscript", "benchmarks/pbmc3k_seurat_github/run_benchmark.R"]
