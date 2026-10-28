# seurust 0.1.0

Fresh baseline release after syncing this fork with upstream **Seurat 5.5.1** and
resetting the seurust version line to **0.1.0** (replacing prior 0.1.x GitHub
release tags).

## Highlights

* **30** native kernels ported to Rust with parity tests against Seurat's C++
  implementations.
* Benchmark tables in the repository `README.md` were refreshed against
  **Seurat 5.5.1**. Per-kernel medians can be reproduced with
  `scripts/bench-all-kernels.R` against an installed copy of Seurat.
* Workflow-level timings (PBMC tutorials and a standard analysis script) used
  the same Seurat 5.5.1 and seurust 0.1.0 pairing. Those helper scripts are
  not part of the installed package.

## Notes

* crates.io may still list older 0.1.x crate versions that were published before
  this reset; they cannot be overwritten. New publishes use this **0.1.0**
  crate version when the registry allows it.
