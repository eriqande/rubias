## Test environments
* local OS X install, R 3.5.1
* ubuntu 14.04 (on travis-ci), R 3.5.2
* win-builder (devel and release)

## R CMD check results

MAC, local: 0 errors | 0 warnings | 0 notes
LINUX, travis-ci: 0 errors | 0 warnings | 1 notes
  The 1 note = libs > 1 Mb due to Rcpp apparently?
WINDOWS, win-builder: 0 errors | 0 warnings | 0 notes

## Downstream dependencies

Currently no known reverse dependencies

## User Notices

Added a number of features (see NEWS) and updated a few things for compatibility with
dplyr 0.8.0.


