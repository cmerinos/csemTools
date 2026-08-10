## Resubmission
This is a new major release (0.2.0) with a new function `curvesCSEMcomp()` for comparing CSEM curves.

## Test environments
* local Windows 11, R 4.5.0
* Ubuntu 22.04 (on GitHub Actions), R 4.4.1, R-devel
* macOS (on GitHub Actions), R 4.4.1

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Reverse dependencies
There are no reverse dependencies.

## Notes for CRAN
This package uses `patchwork` only for optional plotting; it is declared in `Suggests` (moved to `Imports` for CRAN check but still optional). All examples are wrapped in `\donttest{}` and are conditional on required packages being installed.
