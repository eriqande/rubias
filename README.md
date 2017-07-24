rubias --- a package for tidy genetic stock identification (GSI)
================================================================
================
24 July, 2017

-   [Input Data](#input-data)
    -   [An example reference data file](#an-example-reference-data-file)
    -   [An example mixture data file](#an-example-mixture-data-file)
-   [Performing a Genetic Mixture Analysis](#performing-a-genetic-mixture-analysis)
    -   [Aggregating collections into reporting units](#aggregating-collections-into-reporting-units)
    -   [Creating posterior density curves from the traces](#creating-posterior-density-curves-from-the-traces)
    -   [Assessing whether individuals are not from any of the reference populations](#assessing-whether-individuals-are-not-from-any-of-the-reference-populations)
-   [Assessment of Genetic References](#assessment-of-genetic-references)
    -   [Self-assigning fish from the reference](#self-assigning-fish-from-the-reference)
    -   [Simulated mixtures using a leave-one-out type of approach](#simulated-mixtures-using-a-leave-one-out-type-of-approach)
    -   [Specifying mixture proportions in `assess_reference_loo()`](#specifying-mixture-proportions-in-assess_reference_loo)
        -   ["sub-specifying" collection proportions or dirichlet parameters](#sub-specifying-collection-proportions-or-dirichlet-parameters)
        -   ["100% Simulations"](#simulations)
        -   [Do it again with 100% collections](#do-it-again-with-100-collections)
    -   [Bootstrap-Corrected Reporting Unit Proportions](#bootstrap-corrected-reporting-unit-proportions)
-   [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is an R package for performing genetic stock identification (GSI) and associated tasks. Additionally, it includes a method designed to diagnose and correct a bias recently documented in genetic stock identification. The bias occurs when mixture proportion estimates are desired for groups of populations (reporting units) and the number of populations within each reporting unit are uneven.

In order to run C++ implementations of MCMC, rubias requires the package Rcpp, which in turn requires an Rtools installation. After cloning into the repository with the above dependencies installed, build & reload the package to view further documentation.

The script "/R-main/coalescent\_sim" was used to generate coalescent simulations for bias correction validation. This is unnecessary for testing the applicability of our methods to any particular dataset, which should be done using `assess_reference_loo()` and `assess_pb_bias_correction()`. `coalescent_sim()` creates simulated populations using the `ms` coalescent simulation program, available from the Hudson lab at UChicago, and the `GSImulator` and `ms2geno` packages, available at <https://github.com/eriqande>, and so requires more dependencies than the rest of the package.

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
    -   `repunit`: the reporting unit that an individual/collection belongs to. This is required if sample\_type is `reference`. And if sample\_type is `mixture` then repunit must be `NA`.
        This must be a character vector. Not a factor. The idea of a "reporting unit" is well-known amongst people doing genetic stock identfication of salmon, but might not be familiar elsewhere. Briefly, a reporting unit is a group of populations (which we call "collections") that are typically closely related genetically, and which will likely be aggregrated in the results of the GSI exercise.
    -   `collection`: for reference samples, the name of the population that the individual is from. For mixture samples, this is the name of the particular sample (i.e. stratum or port that is to be treated together in space and time.). This must be a character, not a factor.
    -   `indiv` a character vector with the ID of the fish. These must be unique.
-   When we started developing `rubias`, we intended to allow both the `repunit` and the `collection` columns to be either character vectors or factors. Having them as factors might be desirable if, for example, a certain sort order of the collections or repunits was desired. *However* at some point it became clear to Eric that, given our approach to converting all the data to a C++ data structure of integers, for rapid analyis, we would be exposing ourselves to greater opportunities for bugginess by allowing `repunit` and `collection` to be factors. Accordingly, they **must** be character vectors. If they are not, `rubias` will throw an error. **Note**: if you do have a specific sort order for your collections or repunits, you can always change them into factors after analysis with `rubias`. Additionally, you can keep extra columns in your original data frame (for example `repunit_f` or `collection_f`) in which the repunits or the collections are stored as factors. See, for example the data file `alewife`. Or you can just keep a character vector that has the sort order you would like, so as to use it when changing things to factors after `rubias` analysis. (See, for instance, `chinook_repunit_levels`.)
-   The file can have any number of other meta data columns; however, *they must all occur in the data frame **before** the columns of genetic data*.
-   When you pass a data frame into any of these functions, you have to tell it which column the genetic data starts in, and it is assumed that all the columns after that one contain genetic data.
-   If you are doing a mixture analyis, the data frame of mixture fish and of the reference fish must have the same column structure, i.e., they must have exactly the same number of columns with exactly the same column names, in the same order and of the same type.

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

This is done with the `infer_mixture` function. In the example data `chinook_mix` our data consist of fish caught in three different fisheries, `rec1`, `rec2`, and `rec3` as denoted in the collection column. Each of those collections is treated as a separate sample, getting its own mixing proportion estimate. This is how it is run with the default options:

``` r
mix_est <- infer_mixture(reference = chinook, 
                         mixture = chinook_mix, 
                         gen_start_col = 5)
#> Collating data; compiling reference allele frequencies, etc.   time: 7.83 seconds
#> Computing reference locus specific means and variances for computing mixture z-scores   time: 0.72 seconds
#> Working on mixture collection: rec2 with 772 individuals
#>   calculating log-likelihoods of the mixture individuals.   time: 0.29 seconds
#>   performing 100 burn-in and 2000 more sweeps of method "MCMC"   time: 2.03 seconds
#>   tidying output into a tibble.   time: 0.23 seconds
#> Working on mixture collection: rec1 with 743 individuals
#>   calculating log-likelihoods of the mixture individuals.   time: 0.52 seconds
#>   performing 100 burn-in and 2000 more sweeps of method "MCMC"   time: 2.67 seconds
#>   tidying output into a tibble.   time: 0.48 seconds
#> Working on mixture collection: rec3 with 741 individuals
#>   calculating log-likelihoods of the mixture individuals.   time: 0.64 seconds
#>   performing 100 burn-in and 2000 more sweeps of method "MCMC"   time: 2.69 seconds
#>   tidying output into a tibble.   time: 0.23 seconds
```

The result comes back as a list of four tidy data frames:

1.  `mixing_proportions`: the mixing proportions. The column `pi` holds the estimated mixing proportion for each collection.
2.  `indiv_posteriors`: this holds, for each individual, the posterior means of group membership in each collection. Column `PofZ` holds those values. It also includes `n_non_miss_loci` and `n_miss_loci` which are the number of observed loci and the number of missing loci at the individual.
3.  `mix_prop_traces:` MCMC traces of the mixing proportions for each collection.
4.  `bootstrapped_proportions`: This is NULL in the above example, but if we had chosen `method = "PB"` then this would be a tibble of bootstrap-corrected reporting unit mixing proportions.

These data frames look like this:

``` r
lapply(mix_est, head)
#> $mixing_proportions
#> # A tibble: 6 x 4
#>   mixture_collection         repunit           collection           pi
#>                <chr>           <chr>                <chr>        <dbl>
#> 1               rec2 CentralValleyfa         Feather_H_sp 8.008629e-02
#> 2               rec2 CentralValleysp          Butte_Cr_Sp 5.161367e-05
#> 3               rec2 CentralValleysp           Mill_Cr_sp 7.839805e-05
#> 4               rec2 CentralValleysp           Deer_Cr_sp 4.637951e-05
#> 5               rec2 CentralValleysp UpperSacramento_R_sp 4.752110e-04
#> 6               rec2 CentralValleyfa         Feather_H_fa 1.510166e-01
#> 
#> $indiv_posteriors
#> # A tibble: 6 x 10
#>   mixture_collection   indiv         repunit     collection         PofZ
#>                <chr>   <chr>           <chr>          <chr>        <dbl>
#> 1               rec2 T124711 CentralValleyfa   Feather_H_sp 1.843095e-28
#> 2               rec2 T124711 CentralValleyfa   Feather_H_fa 9.998623e-28
#> 3               rec2 T124711 CentralValleyfa    Butte_Cr_fa 1.586087e-24
#> 4               rec2 T124711 CentralValleyfa     Mill_Cr_fa 8.952159e-30
#> 5               rec2 T124711 CentralValleyfa     Deer_Cr_fa 9.876044e-29
#> 6               rec2 T124711 CentralValleyfa Mokelumne_R_fa 1.746479e-27
#> # ... with 5 more variables: log_likelihood <dbl>, z_score <dbl>,
#> #   n_non_miss_loci <dbl>, n_miss_loci <dbl>, missing_loci <list>
#> 
#> $mix_prop_traces
#> # A tibble: 6 x 5
#>   mixture_collection sweep         repunit           collection         pi
#>                <chr> <int>           <chr>                <chr>      <dbl>
#> 1               rec2     0 CentralValleyfa         Feather_H_sp 0.01449275
#> 2               rec2     0 CentralValleysp          Butte_Cr_Sp 0.01449275
#> 3               rec2     0 CentralValleysp           Mill_Cr_sp 0.01449275
#> 4               rec2     0 CentralValleysp           Deer_Cr_sp 0.01449275
#> 5               rec2     0 CentralValleysp UpperSacramento_R_sp 0.01449275
#> 6               rec2     0 CentralValleyfa         Feather_H_fa 0.01449275
#> 
#> $bootstrapped_proportions
#> # A tibble: 0 x 1
#> # ... with 1 variables: mixture_collection <chr>
```

Aggregating collections into reporting units
--------------------------------------------

This is a simple operation in the [tidyverse](http://tidyverse.org/):

``` r
# for mixing proportions
rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit

# for individuals posteriors
rep_indiv_ests <- mix_est$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))
```

Creating posterior density curves from the traces
-------------------------------------------------

The full MCMC output for the mixing proportions is available by default in the field `$mix_prop_traces`. This can be used to obtain an estimate of the posterior density of the mixing proportions.

Here we plot kernel density estimates for the 6 most abundant repunits from the `rec1` fishery:

``` r
# find the top 6 most abundant:
top6 <- rep_mix_ests %>%
  filter(mixture_collection == "rec1") %>% 
  arrange(desc(repprop)) %>%
  slice(1:6)

# check how many MCMC sweeps were done:
nsweeps <- max(mix_est$mix_prop_traces$sweep)

# keep only rec1, then discard the first 200 sweeps as burn-in,
# and then aggregate over reporting units
# and then keep only the top6 from above
trace_subset <- mix_est$mix_prop_traces %>%
  filter(mixture_collection == "rec1", sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top6$repunit)


# now we can plot those:
ggplot(trace_subset, aes(x = repprop, colour = repunit)) +
  geom_density()
```

![](readme-figs/plot-6-1.png)

Assessing whether individuals are not from any of the reference populations
---------------------------------------------------------------------------

Sometimes totally unexpected things happen. One situation we say in the California Chinook fishery was samples coming to us that were actually coho salmon. Before we included coho salmon in the reference sample, these coho always assigned quite strongly to Alaska populations of Chinook, even though they don't really look like Chinook at all.

In this case, it is useful to look at the raw log-likelihood values computed for the individual, rather than the scaled posterior probabilities. Because aberrantly low values of the genotype log-likelihood can indicate that there is something wrong. However, the raw likelihood that you get will depend on the number of missing loci, etc. `rubias` deals with this by computing a *z-score* for each fish. The Z-score is the Z statistic obtained from the fish's log-likelihood (by subtracting from it the expected log-likelihood and dividing by the expected standard deviation). `rubias`'s implementation of the z-score accounts for the pattern of missing data, but it does this without all the simulation that \[`gsi_sim`\])(<https://github.com/eriqande/gsi_sim>) does. This makes it much, much, faster.

Here, we will look at the z-score computed for each fish to the population with the highest posterior. (It is worth noting that you would **never** want to use the z-score to assign fish to different populations---it is only there to decide whether it looks like it reasonably could have come from the population that it was assigned to.)

``` r
# get the maximum-a-posteriori population for each individual
map_rows <- mix_est$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

# then plot them and compare them to bunch of simulated normal rvs
normo <- tibble(z_score = rnorm(1e05))
ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")
```

![](readme-figs/plot-zs-1.png)

The normal density is in black and our observed z\_score is in blue. They fit reasonably well, but it is not perfect. But the z\_score statistic is still a useful check when individuals are wildly aberrant.

Assessment of Genetic References
================================

Self-assigning fish from the reference
--------------------------------------

``` r
sa_chinook <- self_assign(reference = chinook, gen_start_col = 5)
#> Summary Statistics:
#> 
#> 7301 Individuals in Sample
#> 
#> 91 Loci: AldB1.122, AldoB4.183, OTNAML12_1.SNP1, Ots_100884.287, Ots_101119.381, Ots_101704.143, Ots_102213.210, Ots_102414.395, Ots_102420.494, Ots_102457.132, Ots_102801.308, Ots_102867.609, Ots_103041.52, Ots_104063.132, Ots_104569.86, Ots_105105.613, Ots_105132.200, Ots_105401.325, Ots_105407.117, Ots_106499.70, Ots_106747.239, Ots_107074.284, Ots_107285.93, Ots_107806.821, Ots_108007.208, Ots_108390.329, Ots_108735.302, Ots_109693.392, Ots_110064.383, Ots_110201.363, Ots_110495.380, Ots_110551.64, Ots_111312.435, Ots_111666.408, Ots_111681.657, Ots_112301.43, Ots_112419.131, Ots_112820.284, Ots_112876.371, Ots_113242.216, Ots_113457.40, Ots_117043.255, Ots_117242.136, Ots_117432.409, Ots_118175.479, Ots_118205.61, Ots_118938.325, Ots_122414.56, Ots_123048.521, Ots_123921.111, Ots_124774.477, Ots_127236.62, Ots_128302.57, Ots_128693.461, Ots_128757.61, Ots_129144.472, Ots_129170.683, Ots_129458.451, Ots_130720.99, Ots_131460.584, Ots_131906.141, Ots_94857.232, Ots_96222.525, Ots_96500.180, Ots_97077.179, Ots_99550.204, Ots_ARNT.195, Ots_AsnRS.60, Ots_aspat.196, Ots_CD59.2, Ots_CD63, Ots_EP.529, Ots_GDH.81x, Ots_HSP90B.385, Ots_MHC1, Ots_mybp.85, Ots_myoD.364, Ots_Ots311.101x, Ots_PGK.54, Ots_Prl2, Ots_RFC2.558, Ots_SClkF2R2.135, Ots_SWS1op.182, Ots_TAPBP, Ots_u07.07.161, Ots_u07.49.290, Ots_u4.92, OTSBMP.2.SNP1, OTSTF1.SNP1, S71.336, unk_526
#> 
#> 39 Reporting Units: CentralValleyfa, CentralValleysp, CentralValleywi, CaliforniaCoast, KlamathR, NCaliforniaSOregonCoast, RogueR, MidOregonCoast, NOregonCoast, WillametteR, DeschutesRfa, LColumbiaRfa, LColumbiaRsp, MidColumbiaRtule, UColumbiaRsufa, MidandUpperColumbiaRsp, SnakeRfa, SnakeRspsu, NPugetSound, WashingtonCoast, SPugetSound, LFraserR, LThompsonR, EVancouverIs, WVancouverIs, MSkeenaR, MidSkeenaR, LSkeenaR, SSEAlaska, NGulfCoastAlsekR, NGulfCoastKarlukR, TakuR, NSEAlaskaChilkatR, NGulfCoastSitukR, CopperR, SusitnaR, LKuskokwimBristolBay, MidYukon, CohoSp
#> 
#> 69 Collections: Feather_H_sp, Butte_Cr_Sp, Mill_Cr_sp, Deer_Cr_sp, UpperSacramento_R_sp, Feather_H_fa, Butte_Cr_fa, Mill_Cr_fa, Deer_Cr_fa, Mokelumne_R_fa, Battle_Cr, Sacramento_R_lf, Sacramento_H, Eel_R, Russian_R, Klamath_IGH_fa, Trinity_H_sp, Smith_R, Chetco_R, Cole_Rivers_H, Applegate_Cr, Coquille_R, Umpqua_sp, Nestucca_H, Siuslaw_R, Alsea_R, Nehalem_R, Siletz_R, N_Santiam_H, McKenzie_H, L_Deschutes_R, Cowlitz_H_fa, Cowlitz_H_sp, Kalama_H_sp, Spring_Cr_H, Hanford_Reach, PriestRapids_H, Wells_H, Wenatchee_R, CleElum, Lyons_Ferry_H, Rapid_R_H, McCall_H, Kendall_H_sp, Forks_Cr_H, Soos_H, Marblemount_H_sp, QuinaltLake_f, Harris_R, Birkenhead_H, Spius_H, Big_Qual_H, Robertson_H, Morice_R, Kitwanga_R, L_Kalum_R, LPW_Unuk_R, Goat_Cr, Karluk_R, LittleTatsamenie, Tahini_R, Situk_R, Sinona_Ck, Montana_Ck, George_R, Kanektok_R, Togiak_R, Kantishna_R, California_Coho
#> 
#> 4.1752522234648% of allelic data identified as missing
```

Now, you can look at the self assignment results:

``` r
head(sa_chinook, n = 100)
#> # A tibble: 100 x 11
#>                indiv   collection         repunit  inferred_collection
#>                <chr>        <chr>           <chr>                <chr>
#>  1 Feather_H_sp:0001 Feather_H_sp CentralValleyfa         Feather_H_sp
#>  2 Feather_H_sp:0001 Feather_H_sp CentralValleyfa         Feather_H_fa
#>  3 Feather_H_sp:0001 Feather_H_sp CentralValleyfa          Butte_Cr_fa
#>  4 Feather_H_sp:0001 Feather_H_sp CentralValleyfa           Mill_Cr_sp
#>  5 Feather_H_sp:0001 Feather_H_sp CentralValleyfa           Mill_Cr_fa
#>  6 Feather_H_sp:0001 Feather_H_sp CentralValleyfa UpperSacramento_R_sp
#>  7 Feather_H_sp:0001 Feather_H_sp CentralValleyfa           Deer_Cr_sp
#>  8 Feather_H_sp:0001 Feather_H_sp CentralValleyfa          Butte_Cr_Sp
#>  9 Feather_H_sp:0001 Feather_H_sp CentralValleyfa            Battle_Cr
#> 10 Feather_H_sp:0001 Feather_H_sp CentralValleyfa       Mokelumne_R_fa
#> # ... with 90 more rows, and 7 more variables: inferred_repunit <chr>,
#> #   scaled_likelihood <dbl>, log_likelihood <dbl>, z_score <dbl>,
#> #   n_non_miss_loci <dbl>, n_miss_loci <dbl>, missing_loci <list>
```

These can be summarized by repunit as was done above:

``` r
sa_to_repu <- sa_chinook %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))

head(sa_to_repu, n = 200)
#> Source: local data frame [200 x 5]
#> Groups: indiv, collection, repunit [6]
#> 
#> # A tibble: 200 x 5
#>           indiv collection      repunit inferred_repunit repu_scaled_like
#>           <chr>      <chr>        <chr>            <chr>            <dbl>
#>  1 Alsea_R:0001    Alsea_R NOregonCoast  CaliforniaCoast     3.721739e-08
#>  2 Alsea_R:0001    Alsea_R NOregonCoast  CentralValleyfa     1.544394e-14
#>  3 Alsea_R:0001    Alsea_R NOregonCoast  CentralValleysp     8.121327e-15
#>  4 Alsea_R:0001    Alsea_R NOregonCoast  CentralValleywi     1.219934e-23
#>  5 Alsea_R:0001    Alsea_R NOregonCoast           CohoSp     2.092618e-52
#>  6 Alsea_R:0001    Alsea_R NOregonCoast          CopperR     3.077340e-20
#>  7 Alsea_R:0001    Alsea_R NOregonCoast     DeschutesRfa     3.807945e-10
#>  8 Alsea_R:0001    Alsea_R NOregonCoast     EVancouverIs     1.016595e-08
#>  9 Alsea_R:0001    Alsea_R NOregonCoast         KlamathR     1.108037e-11
#> 10 Alsea_R:0001    Alsea_R NOregonCoast     LColumbiaRfa     8.518498e-08
#> # ... with 190 more rows
```

Simulated mixtures using a leave-one-out type of approach
---------------------------------------------------------

If you want to know how much accuracy you can expect given a set of genetic markers and a grouping of populations (`collection`s) into reporting units (`repunit`s), there are two different functions you might use:

1.  `assess_reference_loo()`: This function carries out simulation of mixtures using the leave-one-out approach of Anderson, Waples, and Kalinowski (2008).
2.  `assess_reference_mc()`: This functions breaks the reference data set into different subsets, one of which is used as the reference data set and the other the mixture. It is difficult to simulate very large mixture samples using this method, because it is constrained by the number of fish in the reference data set.
    Additionally, there are constraints on the mixing proportions that can be simulated because of variation in the number of fish from each collection in the reference.

Both of the functions take two required arguments: 1) a data frame of reference genetic data, and 2) the number of the column in which the genetic data start.

Here we use the `chinook` data to simulate 50 mixture samples of size 200 fish using the default values (Dirichlet parameters of 1.5 for each reporting unit, and Dirichlet parameters of 1.5 for each collection within a reporting unit...)

``` r
chin_sims <- assess_reference_loo(reference = chinook, 
                     gen_start_col = 5, 
                     reps = 50, 
                     mixsize = 200)
```

Here is what the output looks like:

``` r
chin_sims
#> # A tibble: 3,450 x 9
#>    repunit_scenario collection_scenario  iter         repunit
#>               <chr>               <chr> <int>           <chr>
#>  1                1                   1     1 CentralValleyfa
#>  2                1                   1     1 CentralValleysp
#>  3                1                   1     1 CentralValleysp
#>  4                1                   1     1 CentralValleysp
#>  5                1                   1     1 CentralValleysp
#>  6                1                   1     1 CentralValleyfa
#>  7                1                   1     1 CentralValleyfa
#>  8                1                   1     1 CentralValleyfa
#>  9                1                   1     1 CentralValleyfa
#> 10                1                   1     1 CentralValleyfa
#> # ... with 3,440 more rows, and 5 more variables: collection <chr>,
#> #   true_pi <dbl>, n <dbl>, post_mean_pi <dbl>, mle_pi <dbl>
```

Specifying mixture proportions in `assess_reference_loo()`
----------------------------------------------------------

By default, each iteration, the proportions of fish from each reporting unit is simulated from a Dirichlet distribution with parameter (1.5,...,1.5). And, within each reporting unit the mixing proportions from different collections are drawn from a Dirichlet distribution with parameter (1.5,...,1.5).

The value of 1.5 for the Dirichlet parameter for reporting units can be changed using the `alpha_repunit`. The Dirichlet parameter for collections can be set using the `alpha_collection` parameter.

Sometimes, however, more control over the composition of the simulated mixtures is desired. This is achieved by passing a two-column *data.frame* to either `alpha_repunit` or `alpha_collection` (or both). If you are passing the data.frame in for `alpha_repunit`, the first column must be named `repunit` and it must contain characters specifying reporting units. In the data.frame for `alpha_collection` the first column must be named `collection` and must hold strings specifying different collections. It is an error if a repunit or collection is specified that does not exist in the reference. However, you do not need to specify a value for every reporting unit or collection. (If they are absent, the value is assumed to be zero.)

The second column of the data frame must be one of `count`, `ppn` or `dirichlet`. These specify, respectively,

1.  the exact count of individuals to be simulated from each repunit (or collection);
2.  the proportion of individuals from each repunit (or collection); or
3.  the parameters of a Dirichlet distribution from which the proportion of individuals should be simulated. These `ppn` values will be normalized to sum to one if they do not. As such, they can be regarded as weights.

Let's say that we want to simulate data that roughly have proportions like what we saw in the Chinook `rec1` fishery. We have those estimates in the variable `top6`:

``` r
top6
#> Source: local data frame [6 x 3]
#> Groups: mixture_collection [1]
#> 
#> # A tibble: 6 x 3
#>   mixture_collection                 repunit    repprop
#>                <chr>                   <chr>      <dbl>
#> 1               rec1         CentralValleyfa 0.82087357
#> 2               rec1                KlamathR 0.06712288
#> 3               rec1                  RogueR 0.06118562
#> 4               rec1         CaliforniaCoast 0.02983790
#> 5               rec1 NCaliforniaSOregonCoast 0.00915214
#> 6               rec1          UColumbiaRsufa 0.00422250
```

We could, if we put those `repprop` values into a `ppn` column, simulate mixtures with exactly those proportions. Or if we wanted to simulate exact numbers of fish in a sample of 345 fish, we could get those values like this:

``` r
round(top6$repprop * 350)
#> [1] 287  23  21  10   3   1
```

and then put them in a `cnts` column.

However, in this case, we want to simulate mixtures that look similar to the one we estimated, but have some variation. For that we will want to supply Dirichlet random variable parmaters in a column named `dirichlet`. If we make the values proportional to the mixing proportions, then, on average that is what they will be. If the values are large, then there will be little variation between simulated mixtures. And if the the values are small there will be lots of variation.
We'll scale them so that they sum to 10---that should give some variation, but not too much. Accordingly the tibble that we pass in as the `alpha_repunit` parameter, which describes the variation in reporting unit proportions we would like to simulate would look like this:

``` r
arep <- top6 %>%
  ungroup() %>%
  mutate(dirichlet = 10 * repprop) %>%
  select(repunit, dirichlet)

arep
#> # A tibble: 6 x 2
#>                   repunit dirichlet
#>                     <chr>     <dbl>
#> 1         CentralValleyfa 8.2087357
#> 2                KlamathR 0.6712288
#> 3                  RogueR 0.6118562
#> 4         CaliforniaCoast 0.2983790
#> 5 NCaliforniaSOregonCoast 0.0915214
#> 6          UColumbiaRsufa 0.0422250
```

Let's do some simulations with those repunit parameters. By default, if we don't specify anything extra for the *collections*, they get dirichlet parameters of 1.5.

``` r
chin_sims_repu_top6 <- assess_reference_loo(reference = chinook, 
                     gen_start_col = 5, 
                     reps = 50, 
                     mixsize = 200,
                     alpha_repunit = arep)
```

Now, we can summarise it by reporting unit and then plot it for the values we are interested in:

``` r
# now, call those repunits that we did not specify in arep "OTHER"
# and then sum up over reporting units
tmp <- chin_sims_repu_top6 %>%
  mutate(repunit = ifelse(repunit %in% arep$repunit, repunit, "OTHER")) %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))
  
```

``` r
# then plot them
ggplot(tmp, aes(x = true_repprop, y = reprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)
```

![](readme-figs/plot-top6-1.png)

Or plot comparing to their "n" value, which is the actual number of fish from each reporting unit in the sample.

``` r
ggplot(tmp, aes(x = repu_n_prop, y = reprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)
```

![](readme-figs/plot-top6-n-1.png)

### "sub-specifying" collection proportions or dirichlet parameters

And let's see how things work if you want some control over the collection proportions in the simulations, too:

``` r
# here for an example, just take the top 20 collections with highest estimated proportions
top20coll_ppn <- mix_est$mixing_proportions %>% 
  filter(mixture_collection == "rec1") %>% 
  arrange(desc(pi)) %>% 
  slice(1:20) %>%
  select(-mixture_collection, -repunit) %>%
  rename(ppn = pi)

top20coll_ppn
#> # A tibble: 20 x 2
#>              collection          ppn
#>                   <chr>        <dbl>
#>  1          Butte_Cr_fa 0.4329420148
#>  2         Feather_H_fa 0.1695781927
#>  3            Battle_Cr 0.1177674882
#>  4         Feather_H_sp 0.0791897463
#>  5       Klamath_IGH_fa 0.0609515373
#>  6         Applegate_Cr 0.0527934140
#>  7                Eel_R 0.0186579776
#>  8       Mokelumne_R_fa 0.0124927306
#>  9            Russian_R 0.0111799183
#> 10              Smith_R 0.0091177120
#> 11        Cole_Rivers_H 0.0083922015
#> 12           Mill_Cr_fa 0.0072397016
#> 13         Trinity_H_sp 0.0061713425
#> 14        Hanford_Reach 0.0041303548
#> 15            Umpqua_sp 0.0020867524
#> 16          Spring_Cr_H 0.0017781357
#> 17        Lyons_Ferry_H 0.0017307058
#> 18      Sacramento_R_lf 0.0011938962
#> 19           Deer_Cr_fa 0.0004697988
#> 20 UpperSacramento_R_sp 0.0003710051
```

NEEDS FURTHER EXPLANATION.

### "100% Simulations"

In the fisheries world, "100% simulations" have been a staple. In these simulations, mixtures are simulated in which 100% of the individuals are from one collection (or reporting unit, I suppose). Eric has never been a big fan of these since they don't necessarily tell you how you might do inferring actual mixtures that you might encounter. Nonetheless, since they have been such a mainstay in the field, it is worthwile showing how to do 100% simulations using `rubias`. Furthermore, when people asked for this feature (thanks Brendan!) it made it clear that Eric had to provide a way to simulate multiple different scenarios without re-processing the reference data set each time. So this is what I came up with: the way we do it is to pass a *list* of scenarios to the `alpha_repunit` or `alpha_collection` option in `assess_reference_loo()`. These can be named lists, if desired. So, for example, let's do 100% simulations for each of the repunits in `arep`:

``` r
arep$repunit
#> [1] "CentralValleyfa"         "KlamathR"               
#> [3] "RogueR"                  "CaliforniaCoast"        
#> [5] "NCaliforniaSOregonCoast" "UColumbiaRsufa"
```

We will let the collections within them just be drawn from a dirichlet distribution with parameter 10 (so, pretty close to equal proportions).

So, to do this, we make a list of data frames with the proportions. We'll give it some names too:

``` r
six_hundy_scenarios <- lapply(arep$repunit, function(x) tibble(repunit = x, ppn = 1.0))
names(six_hundy_scenarios) <- paste("All", arep$repunit, sep = "-")
```

Then, we use it, producing only 5 replicates for each scenario:

``` r
repu_hundy_results <- assess_reference_loo(reference = chinook, 
                     gen_start_col = 5, 
                     reps = 5, 
                     mixsize = 50,
                     alpha_repunit = six_hundy_scenarios,
                     alpha_collection = 10)
#> ++++ Starting in on repunit_scenario All-CentralValleyfa with collection scenario 1 ++++
#> Doing LOO simulations rep 1 of 5
#> Doing LOO simulations rep 2 of 5
#> Doing LOO simulations rep 3 of 5
#> Doing LOO simulations rep 4 of 5
#> Doing LOO simulations rep 5 of 5
#> ++++ Starting in on repunit_scenario All-KlamathR with collection scenario 1 ++++
#> Doing LOO simulations rep 1 of 5
#> Doing LOO simulations rep 2 of 5
#> Doing LOO simulations rep 3 of 5
#> Doing LOO simulations rep 4 of 5
#> Doing LOO simulations rep 5 of 5
#> ++++ Starting in on repunit_scenario All-RogueR with collection scenario 1 ++++
#> Doing LOO simulations rep 1 of 5
#> Doing LOO simulations rep 2 of 5
#> Doing LOO simulations rep 3 of 5
#> Doing LOO simulations rep 4 of 5
#> Doing LOO simulations rep 5 of 5
#> ++++ Starting in on repunit_scenario All-CaliforniaCoast with collection scenario 1 ++++
#> Doing LOO simulations rep 1 of 5
#> Doing LOO simulations rep 2 of 5
#> Doing LOO simulations rep 3 of 5
#> Doing LOO simulations rep 4 of 5
#> Doing LOO simulations rep 5 of 5
#> ++++ Starting in on repunit_scenario All-NCaliforniaSOregonCoast with collection scenario 1 ++++
#> Doing LOO simulations rep 1 of 5
#> Doing LOO simulations rep 2 of 5
#> Doing LOO simulations rep 3 of 5
#> Doing LOO simulations rep 4 of 5
#> Doing LOO simulations rep 5 of 5
#> ++++ Starting in on repunit_scenario All-UColumbiaRsufa with collection scenario 1 ++++
#> Doing LOO simulations rep 1 of 5
#> Doing LOO simulations rep 2 of 5
#> Doing LOO simulations rep 3 of 5
#> Doing LOO simulations rep 4 of 5
#> Doing LOO simulations rep 5 of 5
repu_hundy_results
#> # A tibble: 2,070 x 9
#>       repunit_scenario collection_scenario  iter         repunit
#>                  <chr>               <chr> <int>           <chr>
#>  1 All-CentralValleyfa                   1     1 CentralValleyfa
#>  2 All-CentralValleyfa                   1     1 CentralValleysp
#>  3 All-CentralValleyfa                   1     1 CentralValleysp
#>  4 All-CentralValleyfa                   1     1 CentralValleysp
#>  5 All-CentralValleyfa                   1     1 CentralValleysp
#>  6 All-CentralValleyfa                   1     1 CentralValleyfa
#>  7 All-CentralValleyfa                   1     1 CentralValleyfa
#>  8 All-CentralValleyfa                   1     1 CentralValleyfa
#>  9 All-CentralValleyfa                   1     1 CentralValleyfa
#> 10 All-CentralValleyfa                   1     1 CentralValleyfa
#> # ... with 2,060 more rows, and 5 more variables: collection <chr>,
#> #   true_pi <dbl>, n <dbl>, post_mean_pi <dbl>, mle_pi <dbl>
```

### Do it again with 100% collections

Just to make sure that it is clear how to do this with collections (rather than reporting units) as well, lets do 100% simulations for a handful of the collections. Let's just randomly take 5 of them, and do 6 reps for each:

``` r
set.seed(10)
hundy_colls <- sample(unique(chinook$collection), 5)
hundy_colls
#> [1] "Hanford_Reach" "Applegate_Cr"  "N_Santiam_H"   "Soos_H"       
#> [5] "Feather_H_fa"
```

So, now make a list of those with 100% specifications in the tibbles:

``` r
hundy_coll_list <- lapply(hundy_colls, function(x) tibble(collection = x, ppn = 1.0)) %>%
  setNames(paste("100%", hundy_colls, sep = "_"))
```

Then, do it:

``` r
hundy_coll_results <- assess_reference_loo(reference = chinook, 
                     gen_start_col = 5, 
                     reps = 6, 
                     mixsize = 50,
                     alpha_collection = hundy_coll_list)
hundy_coll_results
#> # A tibble: 2,070 x 9
#>    repunit_scenario collection_scenario  iter         repunit
#>               <chr>               <chr> <int>           <chr>
#>  1                1  100%_Hanford_Reach     1 CentralValleyfa
#>  2                1  100%_Hanford_Reach     1 CentralValleysp
#>  3                1  100%_Hanford_Reach     1 CentralValleysp
#>  4                1  100%_Hanford_Reach     1 CentralValleysp
#>  5                1  100%_Hanford_Reach     1 CentralValleysp
#>  6                1  100%_Hanford_Reach     1 CentralValleyfa
#>  7                1  100%_Hanford_Reach     1 CentralValleyfa
#>  8                1  100%_Hanford_Reach     1 CentralValleyfa
#>  9                1  100%_Hanford_Reach     1 CentralValleyfa
#> 10                1  100%_Hanford_Reach     1 CentralValleyfa
#> # ... with 2,060 more rows, and 5 more variables: collection <chr>,
#> #   true_pi <dbl>, n <dbl>, post_mean_pi <dbl>, mle_pi <dbl>
```

Bootstrap-Corrected Reporting Unit Proportions
----------------------------------------------

These are obtained using `method = "PB"` in `infer_mixture()`. When invoked, this will return the regular MCMC results as before, but also will population the `bootstrapped_proportions` field of the output. Doing so takes a little bit longer, computationally, because there is a good deal of simulation involved:

``` r
mix_est_pb <- infer_mixture(reference = chinook, 
                         mixture = chinook_mix, 
                         gen_start_col = 5,
                         method = "PB")
```

And now we can compare the estimates, showing here the 10 most prevalent repunits, in the `rec1` fishery:

``` r
mix_est_pb$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi)) %>%
  left_join(mix_est_pb$bootstrapped_proportions) %>%
  ungroup() %>%
  filter(mixture_collection == "rec1") %>%
  arrange(desc(repprop)) %>%
  slice(1:10)
```

This is a result that we expect: no appreciable difference, because the reporting units are already very well resolved, so we don't expect that the parametric bootstrap procedure would find any benefit in correcting them.

References
==========

Anderson, Eric C, Robin S Waples, and Steven T Kalinowski. 2008. “An Improved Method for Predicting the Accuracy of Genetic Stock Identification.” *Can J Fish Aquat Sci* 65: 1475–86.
