## Test environments

* local OS X install, R 3.6.2
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

MAC, local: 0 errors | 0 warnings | 1 notes
  1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
LINUX, travis-ci: 0 errors | 0 warnings | 2 notes
  1 note = libs > 1 Mb due to Rcpp apparently
  1 note =  GNU extensions in Makefiles (this is for RcppParallel)
WINDOWS OLDRELEASE, 3.5.3, win-builder: 0 errors | 0 warnings | 1 notes
  1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
WINDOWS RELEASE, 3.6.3, win-builder: 0 errors | 0 warnings | 0 notes
WINDOWS DEVEL, 4.0.0 alpha, win-builder: 0 errors | 0 warnings | 0 notes


## Downstream dependencies

Currently no known reverse dependencies

## User Notices

Changes to fix breaking changes with tibble 3.0.0.

