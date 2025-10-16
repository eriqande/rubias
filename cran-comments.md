## Test environments

* local OS X install, R 4.4.2
* win-builder (devel and release)
* Linux and Windows on RHub with `rhub::check_for_cran()`

## R CMD check results



MAC, local: 0 errors | 1 warnings | 2 notes

  - 1 Warning: 
```
checking CRAN incoming feasibility ... [4s/17s] WARNING
  Maintainer: ‘Eric C. Anderson <eriq@rams.colostate.edu>’
  
  Insufficient package version (submitted: 0.3.4, existing: 0.3.4)
  
  New maintainer:
    Eric C. Anderson <eriq@rams.colostate.edu>
  Old maintainer(s):
    Eric C. Anderson <eric.anderson@noaa.gov>
```

I am the same person, I just retired from federal service and so am using a
different email.  I no longer have access to eric.anderson@noaa.gov.

  - 1 note =  checking HTML version of manual ... NOTE. This seems to be HTML validation problems. These
              are issues I have never seen before. I suspect roxygen2 is at fault as the problems occur
              in all Rd files.  I am using roxygen 7.3.3 which is the most up-to-date.  
  - 1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
  
LINUX, Ubuntu Linux 20.04.1 LTS, R-release, GCC: 0 errors | 0 warnings | 4 notes

  - 1 note = installed size is 10.8Mb, because libs is > 9.0 Mb. This is RCpp on Linux issue I suspect
  - 1 note =  GNU extensions in Makefiles (this is for RcppParallel)
  - 1 note = Examples with CPU elapsed time > 5s. But user time is < 2.2 secs.  I suspect
    the check machine is somewhat oversubscribed.
  - 1 note = Skipping checking HTML validation: no command 'tidy' found (problem with check machine?)
  
LINUX, Fedora Linux, R-devel, clang, gfortran:  0 errors | 0 warnings | 3 notes
  - 1 note =  GNU extensions in Makefiles (this is for RcppParallel)
  - 1 note = Examples with CPU elapsed time > 5s. But user time is < 2.2 secs.  I suspect
    the check machine is somewhat oversubscribed.
  - 1 note = Skipping checking HTML validation: no command 'tidy' found (problem with check machine?)
LINUX, Debian Linux, R-devel, GCC ASAN/UBSAN: 
  - There was a PrepError on this check machine.  Suspect that the config isn't quite right...
WINDOWS RHUB, Windows Server 2022, R-devel, 64 bit: 0 errors | 0 warnings | 3 notes
  - 1 note =  GNU make is a SystemRequirements (this is for RcppParallel)
  - 1 note = Found the following files/directories: ''NULL''  (This seems to be an issue with the check process)
  - 1 note = Found the following files/directories: 'lastMiKTeXException' (PC Latex exception of some sort?)
WINDOWS check_win RELEASE, 4.5.1, win-builder: 0 errors | 0 warnings | 1 notes
  - 1 note = hecking CRAN incoming feasibility because my email address has changed
WINDOWS check_win DEVEL, (2022-02-07 r81667 ucrt), win-builder: 0 errors | 0 warnings | 1 notes
  - 1 note = hecking CRAN incoming feasibility because my email address has changed



## Downstream dependencies

Currently no known reverse dependencies

## revdepcheck results

* 0 reverse dependencies

## User Notices

* This release fixes a few NOTES that Kurt Hornik emailed me about.

