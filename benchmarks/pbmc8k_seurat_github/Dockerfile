ARG BASE_IMAGE=rust-seurat-rust:dev
FROM ${BASE_IMAGE}

WORKDIR /workspace
COPY . /workspace

CMD ["Rscript", "benchmarks/pbmc8k_seurat_github/run_benchmark.R"]
