# Changelog

## [1.0.2] - August 23 2022

### Changed

- Fixed HTML5 validation issue requested by CRAN given R 4.2.0

## [1.0.1] - December 12 2021

### Changed

- Clarify roxygen2 documentation, `palette` arguments
- Fixed `l.max` parameter in `smoothSignalOnGraph`, added validation for graph connectivity
- `extendMatrix` doesn't drop dimensions anymore
- Fixed processing of `mc.allow.recursive` for `n.cores=1` in `plapply`
- Export `heatFilter()`

## [1.0.0] - October 7 2021

### Changed

- The package does not rely on OpenMP anymore, all function use C++11 threads
- Better processing of corner cases in `extendMatrix` and `smoothSignalOnGraph`
- Fixed bug in `val2ggcol` with all color values <= 0: it now produces blue palette instead of read-blue
- Fixed `extendMatrix()` so that always subsets on `col.names`

## [0.1.3] - May 4 2021

### Added

- `collapseCellsByType` from Conos

### Changed

- `colSumByFactor` now requires factor input and returns a named matrix
- Parallel cpp functions can now be accessed from other cpp packages by including `sccore_par.hpp`

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
