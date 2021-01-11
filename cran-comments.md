## Test environments

* local OS X install, R 4.0.3
* ubuntu 16.04.02 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

MAC, local: 0 errors | 0 warnings | 1 notes
  1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
LINUX, travis-ci: 0 errors | 0 warnings | 2 notes
  1 note = libs > 1 Mb due to Rcpp apparently
  1 note =  GNU extensions in Makefiles (this is for RcppParallel)
WINDOWS OLDRELEASE, 3.6.3, win-builder: 0 errors | 0 warnings | 0 notes
WINDOWS RELEASE, 4.0.3, win-builder: 0 errors | 0 warnings | 0 notes
WINDOWS DEVEL, (2021-01-09 r79815), win-builder: 0 errors | 0 warnings | 0 notes


## Downstream dependencies

Currently no known reverse dependencies

## User Notices

* Fixed an underflow issue.
* Fixed the format for DOI's in the DESCRIPTION as requested by Kurt Hornik.

