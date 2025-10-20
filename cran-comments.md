## Test environments

* local OS X install, R 4.4.2
* win-builder (devel and release)
* Linux and Windows on RHub with `rhub::check_for_cran()`

## R CMD check results

## MAC, local: 0 errors | 1 warnings | 1 notes

  - 1 Warning: 
```
checking CRAN incoming feasibility ... [4s/17s] WARNING
  Maintainer: ‘Eric C. Anderson <eriq@rams.colostate.edu>’
  
  New maintainer:
    Eric C. Anderson <eriq@rams.colostate.edu>
  Old maintainer(s):
    Eric C. Anderson <eric.anderson@noaa.gov>
```

I am the same person, I just retired from federal service and so am using a
different email.  I no longer have access to eric.anderson@noaa.gov.
 
  - 1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
```

## win-builder

R Under development (unstable) (2025-10-15 r88920 ucrt)
    1 note = email change (as above)
    
R version 4.5.1 (2025-06-13 ucrt)
    1 note = email change (as above)


## RhubV2
  
RhubV2: ATLAS, R-devel (2025-10-13 r88918), Fedora Linux 38 (Container Image)
  - 0 errors | 0 warnings | 0 notes | 1 INFO = GNU extensions in Makefiles (this is for RcppParallel)

RhubV2: Linux, R-devel (2025-10-13 r88918), Ubuntu 13.3.0
  - 0 errors | 0 warnings | 0 notes | 1 INFO = GNU extensions in Makefiles (this is for RcppParallel)

RhubV2: M1-san, R-devel (2025-10-13 r88918), aarch64-apple-darwin20
  - 0 errors | 0 warnings | 0 notes | 1 INFO = GNU extensions in Makefiles (this is for RcppParallel)

RhubV2: macos, R-devel (2025-10-13 r88918), x86_64-apple-darwin20
  - 0 errors | 0 warnings | 0 notes | 1 INFO = GNU extensions in Makefiles (this is for RcppParallel)
  
RhubV2: macos-arm64, R-devel (2025-10-13 r88918), aarch64-apple-darwin20
  - 0 errors | 0 warnings | 0 notes | 1 INFO = GNU extensions in Makefiles (this is for RcppParallel)

RhubV2: windows, R-devel (2025-10-13 r88918), x86_64-w64-mingw32, Windows Server 2022 x64 (build 26100)
  - 0 errors | 0 warnings | 0 notes | 1 INFO = GNU extensions in Makefiles (this is for RcppParallel)



  

## Downstream dependencies

Currently no known reverse dependencies

## revdepcheck results

* 0 reverse dependencies

## User Notices

* This release adds functionality to sample from the posterior predictive distribution
to propagate uncertainty to total catch.

