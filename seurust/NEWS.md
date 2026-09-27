# seurust 0.1.0

Fresh baseline release after syncing this fork with upstream **Seurat 5.5.1** and
resetting the seurust version line to **0.1.0** (replacing prior 0.1.x GitHub
release tags).

## Highlights

* **30** native kernels ported to Rust with parity tests against Seurat's C++
  implementations.
* Benchmark tables in the root `README.md` refreshed against **Seurat 5.5.1**
  from this repository (`scripts/bench-all-kernels.R`).
* Workflow benchmarks under `benchmarks/` target the same Seurat + seurust
  pairing (PBMC 3K/8K tutorials and `seurat-standard-analysis` scripts).

## Notes

* crates.io may still list older 0.1.x crate versions that were published before
  this reset; they cannot be overwritten. New publishes use this **0.1.0**
  crate version when the registry allows it.
