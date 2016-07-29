# rubias --- a package for bias correction in hierarchical GSI

This is an R package designed to diagnose and correct a bias
which appears in genetic stock identification, when mixture proportion
estimates are desired for groups of populations (here called reporting units)
and the number of populations within each reporting unit are uneven.

In order to run C++ implementations of MCMC, rubias requires the package
Rcpp, which in turn requires an Rtools installation. After cloning into the
repository with the above dependencies installed, build & reload the package 
to view further documentation.

The script "/R-main/coalescent_sim" was used to generate coalescent simulations
for bias correction validation. This is unnecessary for testing the applicability
of our methods to any particular dataset, which should be done using the 
Hasselman_simulation_pipeline and/or bias_comparison functions.
coalescent_sim creates simulated populations using the ms coalescent simulation 
program, available from the Hudson lab at UChicago, and the GSImulator
and ms2geno packages, available at https://github.com/eriqande, and so requires
further dependencies than the rest of the package.
