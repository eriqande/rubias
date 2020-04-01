

library(tidyverse)
library(rubias)

#### Chinook Speed Test On Old Mac with 24 cores ####

ref <- read.csv("CI_Chinook_67pops_429loci_6repunits_base.csv", stringsAsFactors = FALSE, colClasses = c("character"))
mix <- read.csv("GeneralS17_mix.csv", stringsAsFactors = FALSE, colClasses = c("character"))

system.time(plain <- infer_mixture(reference = ref, mixture = mix, gen_start_col = 5))


# Collating data; compiling reference allele frequencies, etc.   time: 29.31 seconds
# Computing reference locus specific means and variances for computing mixture z-scores   time: 3.71 seconds
# Working on mixture collection: GeneralS17 with 426 individuals
# calculating log-likelihoods of the mixture individuals.   time: 0.71 seconds
# performing 2000 total sweeps, 100 of which are burn-in and will not be used in computing averages in method "MCMC"   time: 0.65 seconds
# tidying output into a tibble.   time: 0.09 seconds
#
#
# user  system elapsed
# 58.569   3.105  61.852


system.time(br_res <- infer_mixture(reference = ref, mixture = mix, gen_start_col = 5, method = "BR"))

# Collating data; compiling reference allele frequencies, etc.   time: 22.43 seconds
# Computing reference locus specific means and variances for computing mixture z-scores   time: 3.71 seconds
# Working on mixture collection: GeneralS17 with 426 individuals
# calculating log-likelihoods of the mixture individuals.   time: 0.71 seconds
# performing 2000 sweeps of method "BR", 100 sweeps of which are burn-in.   time: 94.90 seconds
# tidying output into a tibble.   time: 0.09 seconds
#
# user   system  elapsed
# 1226.815   43.446  149.198


#### Chinook Speed Test on Linux Box with 12 cores ####

ref <- read.csv("CI_Chinook_67pops_429loci_6repunits_base.csv", stringsAsFactors = FALSE, colClasses = c("character"))
mix <- read.csv("GeneralS17_mix.csv", stringsAsFactors = FALSE, colClasses = c("character"))

system.time(plain <- infer_mixture(reference = ref, mixture = mix, gen_start_col = 5))

# Collating data; compiling reference allele frequencies, etc.   time: 19.50 seconds
# Computing reference locus specific means and variances for computing mixture z-scores   time: 2.46 seconds
# Working on mixture collection: GeneralS17 with 426 individuals
# calculating log-likelihoods of the mixture individuals.   time: 1.03 seconds
# performing 2000 total sweeps, 100 of which are burn-in and will not be used in computing averages in method "MCMC"   time: 0.72 seconds
# tidying output into a tibble.   time: 0.08 seconds
#
# user  system elapsed
# 25.215   1.036  26.415


system.time(br_res <- infer_mixture(reference = ref, mixture = mix, gen_start_col = 5, method = "BR"))

# Collating data; compiling reference allele frequencies, etc.   time: 15.36 seconds
# Computing reference locus specific means and variances for computing mixture z-scores   time: 2.39 seconds
# Working on mixture collection: GeneralS17 with 426 individuals
# calculating log-likelihoods of the mixture individuals.   time: 1.03 seconds
# performing 2000 sweeps of method "BR", 100 sweeps of which are burn-in.   time: 180.86 seconds
# tidying output into a tibble.   time: 0.14 seconds
#
# user   system  elapsed
# 1766.400    2.320  202.533
