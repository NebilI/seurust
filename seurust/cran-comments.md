## Resubmission

This is a resubmission as **0.1.2** with an updated maintainer email address.

Previous 0.1.0 incoming pretests failed with:

* Windows ERROR: missing `configure.win`, so `Makevars.win` was never
  generated and only `entrypoint.c` was linked (undefined Rust symbols).
  Fixed in 0.1.1 by adding `configure.win` / `cleanup.win`.
* Debian WARNING: GNU make extension `.NOTPARALLEL` in
  `src/Makevars` / `src/Makevars.in`. Removed; cleanup remains ordered via
  ordinary Make dependencies (`rust_clean: $(SHLIB)`).
* NOTE: DESCRIPTION wording rewritten to avoid spell-check false positives;
  README no longer links to package-local `CRAN.md` (that file is
  `.Rbuildignore`d).

## Test environments

* GitHub Actions (ubuntu-latest), R release — `build-seurust-cran.yaml`
* Prior Docker CRAN check for the 0.1.1 packaging fixes (Status: 1 NOTE —
  New submission only)

## R CMD check results

There were no ERRORs.

Notes expected for this package:

* Compiled code uses Rust plus a small C++ bridge for ModularityOptimizer.
  `SystemRequirements` lists Cargo and rustc (>= 1.81).
* Source tarball includes offline Rust crate sources in
  `src/rust/vendor.tar.xz` (~2 MB) so CRAN can build offline.

## Downstream dependencies

There are currently no reverse dependencies on CRAN.
