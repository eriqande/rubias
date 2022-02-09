# rubias 0.3.3

## Changes

* Overhauled haploid vs diploid ploidy determination to deal gracefully with mixture samples that
are missing data at all individuals in a locus.
* Using `dplyr::slice_sample()` instead of `dplyr::sample_n()`, as the latter is superseded
by the former, and the latter was also causing an error on the new sample.int() sanity check.



# rubias 0.3.2

## Changes

* Corrected an underflow issues in the self_assign() function with
many loci.
* Edited the DOI address of the paper to be CRAN compliant.

# rubias 0.3.1

## Changes

* Tiny tweak in one function to fix a breaking change in tibble 3.0.0.

# rubias 0.3.0

## Changes

* Added option to `infer_mixture()` for a fully-Bayesian version.
In a fully Bayesian version, fish from within
the mixture that are allocated (on any particular step of the MCMC) to one of the reference
samples have their alleles added to that reference sample.


# rubias 0.2.0

## Changes

* Added support for haploid markers (#14, @krshedd).

* Added support for individuals of known origin (i.e. those identified with great accuracy using 
parentage-based tagging) in the mixtures (#12).

* Allow user to specify the total weight on the symmetrical Dirichlet prior for the mixing
proportions in infer_mixture().

* Enforced the requirement that fish of sample_type == "mixture" must have NA for their repunit.
When things aren't NA, infer_mixture() would throw an error when method == "PB" because there 
were extra factor levels floating around.  In the process I allow for the repunit column to be
either character or logical as setting it to NA will always make it a logical if it is not part
of a data frame with other non-missing character values in it.

* Modified vignette so that we don't have to put tidyverse in the Suggests

* Added support for user specified parameters for the Dirichlet prior on mixing proportions

* Added support for user-specified initial starting values for the pi parameter (the mixing proportions of collections)
in the function infer_mixture().

* Added a simple function, close_matching_samples(), to tabulate pairs of
individuals with only a small number of mismatching
genotypes.  This is useful for identifying accidentally duplicated samples. (#23)

* Made changes to be compatible with dplyr 0.8.0, which no longer discards
empty factor levels.  Mostly this involved filtering 0's after group_by() and
tally()  or count() calls.  It looks like things are all working (passing all tests
and building vignettes all right.)



# rubias 0.1.0

This is the first version submitted to CRAN.
