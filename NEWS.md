
# rubias 0.1.0.9000

## Changes

* Added support for haploid markers (#14, @krshedd).

* Added support for individuals of known origin (i.e. those identified with great accuracy using 
parentage-based tagging) in the mixtures (#12).

* Allow user to specify the total weight on the symmetrical Dirichlet prior for the mixing
proportions in infer_mixture().

* Enforced the requirement that fish of sample_type == "mixture" must have NA for their repunit.
When things aren't NA, infer_mixture() would throw an error when method == "PB" because there 
were extra factor levels floating around.  This 



# rubias 0.1.0

This is the first version submitted to CRAN.
