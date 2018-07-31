
# rubias 0.1.0.9000

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





# rubias 0.1.0

This is the first version submitted to CRAN.
