## Test environments
* local OS X install, R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

## R CMD check results

MAC, local: 0 errors | 0 warnings | 1 notes
  1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
LINUX, travis-ci: 0 errors | 0 warnings | 2 notes
  1 note = libs > 1 Mb due to Rcpp apparently
  1 note =  GNU extensions in Makefiles (this is for RcppParallel)
WINDOWS, win-builder: 0 errors | 0 warnings | 0 notes
  1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
  
## Downstream dependencies

Currently no known reverse dependencies

## User Notices

<<<<<<< HEAD
Added a number of features (see NEWS) and updated a few things for compatibility with
dplyr 0.8.0.
=======
Added a "fully-Bayesian" option to the infer_mixture() function.
>>>>>>> baseline_resampling


