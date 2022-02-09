## Test environments

* local OS X install, R 4.1.1
* win-builder (devel and release)
* Linux on Rhub

## R CMD check results

MAC, local: 0 errors | 0 warnings | 1 notes
  - 1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
LINUX, Ubuntu 20.04: 0 errors | 0 warnings | 3 notes
  - 1 note = checking CRAN incoming feasibility, X-CRAN-Comment: Archived on
    2022-02-06 as check problems were not
    corrected in time.
  - 1 note = libs > 1 Mb due to Rcpp apparently
  - 1 note =  GNU extensions in Makefiles (this is for RcppParallel)
LINUX, Fedora:  0 errors | 0 warnings | 3 notes
  - 1 note = checking CRAN incoming feasibility, X-CRAN-Comment: Archived on
    2022-02-06 as check problems were not
    corrected in time.
  - 1 note = libs > 1 Mb due to Rcpp apparently
  - 1 note =  GNU extensions in Makefiles (this is for RcppParallel)
WINDOWS RELEASE, 4.1.2, win-builder: 0 errors | 0 warnings | 1 notes
  - note = checking incoming cran feasibility: X-CRAN-Comment: Archived on 
    2022-02-06 as check problems were not corrected in time.
WINDOWS DEVEL, (2022-02-07 r81667 ucrt), win-builder: 0 errors | 1 warnings | 0 notes
  - 1 warning = checking incoming cran feasibility: X-CRAN-Comment: Archived on 
    2022-02-06 as check problems were not corrected in time.


## Downstream dependencies

Currently no known reverse dependencies

## User Notices

* This was archived from CRAN on 2022-02-06 because it didn't pass the new sample.int() sanity
check.  The problem occurs in a call to one of dplyr's functions `sample_n()`. That function
has been superseded in the dplyr package by `slice_sample()`.  I replaced `sample_n()` with 
`slice_sample()` and this now passes the checks, including windows devel
(which errored on the old version using `sample_n()`).

