ARG BASE_IMAGE=seurust-rust:dev
FROM ${BASE_IMAGE}

WORKDIR /workspace
COPY . /workspace

CMD ["Rscript", "benchmarks/seurat_standard_analysis_github/run_benchmark.R"]
