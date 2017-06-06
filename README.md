rubias --- a package for bias correction in hierarchical GSI
============================================================
================
05 June, 2017

-   [Input Data](#input-data)
    -   [An example reference data file](#an-example-reference-data-file)
    -   [An example mixture data file](#an-example-mixture-data-file)
-   [Performing a Genetic Mixture Analysis](#performing-a-genetic-mixture-analysis)
-   [Simulation Assessment of Genetic References](#simulation-assessment-of-genetic-references)
    -   [Specifying mixture proportions in `assess_reference_loo()`](#specifying-mixture-proportions-in-assess_reference_loo)
-   [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is an R package for perorming genetic stock identification (GSI) and associated tasks. Additionally, it includes a method designed to diagnose and correct a bias recently documented in genetic stock identification. The bias occurs when mixture proportion estimates are desired for groups of populations (reporting units) and the number of populations within each reporting unit are uneven.

In order to run C++ implementations of MCMC, rubias requires the package Rcpp, which in turn requires an Rtools installation. After cloning into the repository with the above dependencies installed, build & reload the package to view further documentation.

The script "/R-main/coalescent\_sim" was used to generate coalescent simulations for bias correction validation. This is unnecessary for testing the applicability of our methods to any particular dataset, which should be done using `assess_reference_loo()` and `assess_bp_bias_correction()`. `coalescent_sim()` creates simulated populations using the `ms` coalescent simulation program, available from the Hudson lab at UChicago, and the `GSImulator` and `ms2geno` packages, available at <https://github.com/eriqande>, and so requires more dependencies than the rest of the package.

Input Data
==========

The functions for conducting genetic mixture analysis and for doing simulation assessment to predict the accuracy of a set of genetic markers for genetic stock identification require that genetic data be input as a data frame in a specific format:

-   one row per individual
-   each locus is represented by two adjacent columns, one for each allele (this package is only configured for diploids, at the moment). Allelic types can be expressed as any number or character
-   missing data at a locus is expressed with NA values for each gene copy at the locus
-   if one gene copy is missing in an indivividuals from a locus, then both gene copies must be missing at the locus.
-   the name of the locus is taken to be the column name of the *first* column of each pair of locus columns. The header on the second column is ignored.
-   the data frame must have four columns of meta data for each individual:
    -   `sample_type`: a column telling whether the sample is a `reference` sample or a `mixture` sample.
    -   `repunit`: the reporting unit that an individual belongs to. This is required if sample\_type is `reference`. And if sample\_type is `mixture` then repunit must be `NA`. This must be a character vector. Not a factor.
    -   `collection`: for reference samples, the name of the population that the individual is from. For mixture samples, this is the name of the particular sample (i.e. stratum or port that is to be treated together in space and time.). This must be a character, not a factor.
    -   `indiv` a character vector with the ID of the fish. These must be unique.
-   When we started developing `rubias`, we intended to allow both the `repunit` and the `collection` columns to be either character vectors or factors. Having them as factors might be desirable if, for example, a certain sort order of the collections or repunits was desired. *However* at some point it became clear to Eric that, given our approach to converting all the data to a C++ data structure of integers, we would be exposing ourselves to greater opportunities for bugginess by allowing `repunit` and `collection` to be factors. Accordingly, they **must** be character vectors. If they are not, `rubias` will throw an error. **Note**: if you do have a specific sort order for your collections or repunits, you can always change them into factors after analysis with `rubias`. Additionally, you can keep extra columns in your original data frame (for example `repunit_f` or `collection_f`) in which the repunits or the collections are stored as factors. See, for example the data file `alewife`
-   The file can have any number of other meta data columns; however, *they must all occur in the data frame **before** the columns of genetic data*.
-   When you pass a data frame into any of these functions, you have to tell it which column the genetic data starts in, and it is assumed that all the columns after that one contain genetic data.
-   If you are doing a mixture analyis, the data frame of mixture fish and of the reference fish must have the same column structure, i.e., they must have exactly the same number of columns with exactly the same column names, in the same order.

An example reference data file
------------------------------

Here are the meta data columns and the first two loci for eight individuals in the `chinook` reference data set that comes with the package:

``` r
library(rubias)
head(chinook[, 1:8])
#> # A tibble: 6 x 8
#>   sample_type         repunit   collection             indiv Ots_94857.232
#>         <chr>           <chr>        <chr>             <chr>         <int>
#> 1   reference CentralValleyfa Feather_H_sp Feather_H_sp:0001             2
#> 2   reference CentralValleyfa Feather_H_sp Feather_H_sp:0002             2
#> 3   reference CentralValleyfa Feather_H_sp Feather_H_sp:0003             2
#> 4   reference CentralValleyfa Feather_H_sp Feather_H_sp:0004             2
#> 5   reference CentralValleyfa Feather_H_sp Feather_H_sp:0005             2
#> 6   reference CentralValleyfa Feather_H_sp Feather_H_sp:0006             2
#> # ... with 3 more variables: Ots_94857.232.1 <int>, Ots_102213.210 <int>,
#> #   Ots_102213.210.1 <int>
```

An example mixture data file
----------------------------

Here is the same for the mixture data frame that goes along with that reference data set:

``` r
head(chinook_mix[, 1:8])
#> # A tibble: 6 x 8
#>   sample_type repunit collection   indiv Ots_94857.232 Ots_94857.232.1
#>         <chr>   <chr>      <chr>   <chr>         <int>           <int>
#> 1     mixture    <NA>       rec2 T124711             4               2
#> 2     mixture    <NA>       rec2 T124719             4               2
#> 3     mixture    <NA>       rec2 T124727             4               4
#> 4     mixture    <NA>       rec1 T124735             4               4
#> 5     mixture    <NA>       rec1 T124743             2               2
#> 6     mixture    <NA>       rec1 T124759             4               2
#> # ... with 2 more variables: Ots_102213.210 <int>, Ots_102213.210.1 <int>
```

Performing a Genetic Mixture Analysis
=====================================

Not written.

Simulation Assessment of Genetic References
===========================================

If you want to know how much accuracy you can expect given a set of genetic markers and a grouping of populations (`collection`s) into reporting units (`repunit`s), there are two different functions you might use:

1.  `assess_reference_loo()`: This function carries out simulation of mixtures using the leave-one-out approach of Anderson, Waples, and Kalinowski (2008).
2.  `assess_reference_mc()`: This functions breaks the reference data set into different subsets, one of which is used as the reference data set and the other the mixture. It is difficult to simulate very large mixture samples using this method, because it is constrained by the number of fish in the reference data set.
    Additionally, there are constraints on the mixing proportions that can be simulated because of variation in the number of fish from each collection in the reference.

Both of the functions take two required arguments: 1) a data frame of reference genetic data, and 2) the number of the column in which the genetic data start.

Specifying mixture proportions in `assess_reference_loo()`
----------------------------------------------------------

By default, each iteration, the proportions of fish from each reporting unit is simulated from a Dirichlet distribution with parameter (1.5,...,1.5). And, within each reporting unit the mixing proportions from different collections are drawn from a Dirichlet distribution with parameter (1.5,...,1.5).

The value of 1.5 for the Dirichlet parameter for reporting units can be changed using the `alpha_repunit`. The Dirichlet parameter for collections can be set using the `alpha_collection` parameter.

Sometimes, however, more control over the composition of the simulated mixtures is desired. This is achieved by passing a two-column *data.frame* to either `alpha_repunit` or `alpha_collection` (or both). If you are passing the data.frame in for `alpha_repunit`, the first column must be named `repunit` and it must contain characters specifying reporting units. In the data.frame for `alpha_collection` the first column must be named `collection` and must hold strings specifying different collections. It is an error if a repunit or collection is specified that does not exist in the reference. However, you do not need to specify a value for every reporting unit or collection. (If they are absent, the value is assumed to be zero.)

The second column of the data frame must be one of `count`, `ppn` or `dirichlet`. These specify, respectively, 1) the exact count of individuals to be simulated from each repunit (or collection); 2) the proportion of individuals from each repunit (or collection); or 3) the parameters of a Dirichlet distribution from which the proportion of individuals should be simulated. These `ppn` values will be normalized to sum to one if they do not. As such, they can be regarded as weights.

References
==========

Anderson, Eric C, Robin S Waples, and Steven T Kalinowski. 2008. “An Improved Method for Predicting the Accuracy of Genetic Stock Identification.” *Can J Fish Aquat Sci* 65: 1475–86.
