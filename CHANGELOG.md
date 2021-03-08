## Upcoming

### Added

### Changed

- colSumByFactor now requires factor input and returns a named matrix

## [0.1.2] - February 19 2021

### Added

- `val2col` function that translates values into colors
- added `smoothSignalOnGraph` function that re-implements graph filtering from the [pygsp](https://github.com/epfl-lts2/pygsp/) package
- Rcpp functions `runTaskParallelFor`, `runTaskParallel` and class `ThreadProgress` for parallel loops with progress bar using C++11 threads
- roxygen2 revisions

### Changed

- `plot.na` can accept numeric values now. If plot.na passed a numeric value below 0, the NA symbols are plotted below the cells. Otherwise if values >=0, they're plotted above the cells.
- `plapply` now accepts `fail.on.error` parameter that is `FALSE` by default
- Small bug fixes
- `plapply` now uses `pbmcapply`

## [0.1.1] - September 30 2020
