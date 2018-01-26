## Test environments
* local OS X install, R 3.4.3
* ubuntu 14.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results

MAC, local: 0 errors | 0 warnings | 0 notes
LINUX, travis-ci: 0 errors | 0 warnings | 1 notes
  The 1 note = libs > 1 Mb due to Rcpp apparently?
WINDOWS, win-builder: 0 errors | 0 warnings | 1 note
  The 1 note = first submission

## Downstream dependencies

Currently no known reverse dependencies

## User Notices

This is the first submission of this package.  Swetlana Herbrandt kindly provided
comments yesterday when I first submitted this.  Pursuant to her suggestions, 
I have:

1. Used `person()` in the form `c(person(), person())` in DESCRIPTION, 
rather than `as.person(c(...))`

2. Removed all \dontruns (except one), so that all of the functions get tested
when running examples in in R CMD CHECK.  The one that is still a dontrun
is computationally intensive and can't be done in < 5 sec.

3. `write_gsi_sim_reference()` and `write_gsi_sim_mixture()` both require user-specified
paths for writing a file.  There is no default.  These functions now have examples
which write to files in a `tempdir()` directory.

Thank you very much.

