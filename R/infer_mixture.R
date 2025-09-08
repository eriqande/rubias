#' Estimate mixing proportions and origin probabilities from one or several mixtures
#'
#' Takes a mixture and reference dataframe of two-column genetic data, and a
#' desired method of estimation for the population mixture proportions (MCMC, PB, BR).
#' Returns the output of the chosen estimation method
#'
#' "MCMC" estimates mixing proportions and individual posterior
#' probabilities of assignment through Markov-chain Monte Carlo
#' conditional on the reference allele frequencies,
#' while "PB" does the same with a parametric bootstrapping correction,
#' and "BR" runs MCMC sweeps while simulating reference allele frequencies
#' using the genotypes of mixture individuals and allocations from the previous sweep.
#' All methods default to a uniform 1/(# collections or RUs) prior for the mixing proportions.
#'
#' @param reference a dataframe of two-column genetic format data, proceeded by "repunit", "collection",
#' and "indiv" columns. Does not need "sample_type" column, and will be overwritten if provided
#' @param mixture a dataframe of two-column genetic format data. Must have the same structure as
#' \code{reference} dataframe, but "collection" and "repunit" columns are ignored.
#' Does not need "sample_type" column, and will be overwritten if provided
#' @param gen_start_col the first column of genetic data in both data frames
#' @param method a choice between "MCMC", "PB", "BR" methods for estimating mixture proportions
#' @param alle_freq_prior a one-element named list specifying the prior to be used when
#' generating Dirichlet parameters for genotype likelihood calculations. Valid methods include
#' \code{"const"}, \code{"scaled_const"}, and \code{"empirical"}. See \code{?list_diploid_params}
#' for method details.
#' @param pi_prior The prior to be added to the collection allocations, in order to generate pseudo-count
#' Dirichlet parameters for the simulation of new pi vectors in MCMC. Default value of NA leads to the
#' calculation of a symmetrical prior based on \code{pi_prior_sum}. To provide other values to
#' certain collections, you can pass in a data frame with two columns, "collection"
#' listing the relevant collection, and "pi_param" listing the desired prior for that collection.
#' Specific priors may be listed for as few as one collection. The special collection name "DEFAULT_PI"
#' is used to set the prior for all collections not explicitly listed; if no "DEFAULT_PI" is given, it is
#' taken to be 1/(# collections).
#' @param pi_init  The initial value to use for the mixing proportion of collections.  This lets
#' the user start the chain from a specific value of the mixing proportion vector.  If pi_init is NULL
#' (the default) then the mixing proportions are all initialized to be equal.  Otherwise, you pass
#' in a data frame with one column named "collection" and the other named "pi_init".  Every value in the
#' pi_init column must be strictly positive (> 0), and a value must be given for every collection.  If they sum
#' to more than one the values will be normalized to sum to one.
#' @param reps the number of iterations to be performed in MCMC
#' @param burn_in how many reps to discard in the beginning of MCMC when doing the mean calculation.
#' They will still be returned in the traces if desired.
#' @param pb_iter how many bootstrapped data sets to do for bootstrap correction using method PB.  Default
#' is 100.
#' @param prelim_reps for method "BR", the number of reps of conditional MCMC (as in method "MCMC")
#' to perform prior to MCMC with baseline resampling. The posterior mean of mixing proportions
#' from this conditional MCMC is then used as \code{pi_init} in the baseline resampling MCMC.
#' @param prelim_burn_in for method "BR", this sets the number of sweeps out of \code{prelim_reps}
#' that should be discarded as burn in when preparing the posterior means of the mixing
#' proportions to be set as \code{pi_init} in the baseline resampling MCMC.
#' @param sample_int_Pi how many iterations between storing the mixing proportions trace. Default is 1.
#' Can't be 0. Can't be so large that fewer than 10 samples are taken from the burn in and the sweeps.
#' @param sample_theta for method "BR", whether or not the function should store the posterior mean
#' of the updated allele frequences. Default is TRUE
#' @param pi_prior_sum For \code{pi_prior = NA}, the prior on the mixing proportions is set
#' as a Dirichlet vector of length C, with each element being W/C, where W is the pi_prior_sum
#' and C is the number of collections. By default this is 1.  If it is made much smaller than 1, things
#' could start to mix more poorly.
#' @param total_catch_tib  A tibble with two columns: `collection` and `tot_catch`. The entries to `collection`
#' must correspond to different fisheries or mixture collections as listed in the `collections` column in
#' `mixture`.  The `tot_catch` must be a list column which contains an estimate of the total
#' catch from the fishery.  That could be a single number, or it could be a sample from a posterior
#' distribution of the total catch.  Providing a value for this parameter will make the function
#' return a sample from the posterior for stock-specific total catch. Note that if you use this option
#' then every single mixture collection in `mixture` must be represented in `total_catch_tib`.
#' @param variable_prob_is_catch a logical indicating whether to use the `catch_is_prob` column to simulate
#' for each fish whether or not it should be considered part of the catch.
#' @return Tidy data frames in a list with the following components:
#' mixing_proportions: the estimated mixing proportions of the different collections.
#' indiv_posteriors: the posterior probs of fish being from each of the collections.
#' mix_prop_traces: the traces of the mixing proportions.  Useful for computing credible intervals.
#' bootstrapped_proportions: If using method "PB" this returns the bootstrap corrected
#' reporting unit proportions.  If `total_catch_tib` is not NULL, then you will also found
#' three other components in the output list: `stock_specific_total_catch_traces`,
#' `posterior_predictive_remaining_catch_traces` and `allocation_count_traces`.  These
#' give samples from the posterior distribution of total catch.  Note that
#' `stock_specific_total_catch_traces` is the sum of `posterior_predictive_remaining_catch_traces`
#' and `allocation_count_traces`.  All three are returned because it can be instructive
#' to be able to decompose the `stock_specific_total_catch_traces` into its two
#' constituents.
#'
#' @examples
#' mcmc <- infer_mixture(reference = small_chinook_ref,
#'                       mixture = small_chinook_mix,
#'                       gen_start_col = 5,
#'                       method = "MCMC",
#'                       reps  = 200)
#' @export
infer_mixture <- function(reference,
                          mixture,
                          gen_start_col,
                          method = "MCMC",
                          alle_freq_prior = list("const_scaled" = 1),
                          pi_prior = NULL,
                          pi_init = NULL,
                          reps = 2000,
                          burn_in = 100,
                          pb_iter = 100,
                          prelim_reps = NULL,
                          prelim_burn_in = NULL,
                          sample_int_Pi = 1,
                          sample_theta = TRUE,
                          pi_prior_sum = 1,
                          total_catch_tib = NULL,
                          variable_prob_is_catch = FALSE
                        ) {


  # check that prelim_reps and prelim_burn_in are appropriately set
  if(xor(is.null(prelim_reps), is.null(prelim_burn_in))) {
    stop("Both prelim_reps and prelim_burn_in must be set.  Not just one of them.")
  }
  if(!is.null(prelim_reps) && !is.null(prelim_burn_in)) {
    if(prelim_reps <= prelim_burn_in) stop("prelim_reps is total reps before discarding prelim_burn_in reps, and hence prelim_reps must be strictly larger than prelim_burn_in.")
  }

  # check that reference and mixture are OK
  ploidies_ref <- check_refmix(reference, gen_start_col, "reference")
  ploidies_mix <- check_refmix(mixture, gen_start_col, "mixture")


  # check that total_catch_tib is valid, if it exists
  if(!is.null(total_catch_tib)) {
    # ensure collection and tot_catch columns exist
    tcm <- setdiff(c("collection", "tot_catch"), names(total_catch_tib))
    if(length(tcm) > 0) stop("total_catch_tib missing necessary column(s): ", paste(tcm, collapse = ", "))
    # ensure tot_catch is a list-column
    if(class(total_catch_tib$tot_catch) != "list") stop("Column tot_catch of total_catch_tib must be a list-column")
    # ensure that a total catch is listed for every collection in mixture
    tc_missing_colls <- setdiff(unique(mixture$collection), unique(total_catch_tib$collection))
    if(length(tc_missing_colls) > 0) stop("total_catch_tib missing entries for collection(s) ", paste(tc_missing_colls, collapse = ", "))
    # finally, ensure that total_catch_tib does not have additional collections not found in mixture
    novel_coll <- setdiff(unique(total_catch_tib$collection), unique(mixture$collection))
    if(length(novel_coll) > 0) stop("total_catch_tib collection(s) not included in mixture: ", paste(novel_coll, collapse = ", "))
  }

  # If variable_prob_is_catch is TRUE, we have to ensure several things:
  # 1. make sure total_catch_tib is not NULL
  # 2. make sure that the reference and mixture genotypes tibbles have a column
  #    named prob_is_catch, and that every value is between 0 and 1 inclusive
  #    in the mixture version.
  if(variable_prob_is_catch == TRUE) {
    if(is.null(total_catch_tib)) {
      stop("When using variable_prob_is_catch = TRUE, you must also supply a total_catch_tib.")
    }
    if(!("prob_is_catch" %in% names(reference))) {
      stop("When using variable_prob_is_catch = TRUE, there must be a column in reference named prob_is_catch.")
    }
    if(!("prob_is_catch" %in% names(mixture))) {
      stop("When using variable_prob_is_catch = TRUE, there must be a column in mixture named prob_is_catch.")
    }
    if(!all(mixture$prob_is_catch >= 0 & mixture$prob_is_catch <= 1)) {
      stop("When using variable_prob_is_catch = TRUE, all entries in the prob_is_catch column in mixture must be between 0 and 1 inclusive.")
    }
  }


  # now, deal with the problem if ploidy is indeterminate in both ref and mix
  both_indet <- ploidies_ref == 0 & ploidies_mix == 0
  if(any(both_indet == TRUE)) {
    stop(
      "All allelic data missing for the following loci in both the reference and the mixture. Bailing out. ",
      paste(names(both_indet)[both_indet], collapse = ", ")
    )
  }

  # also give a warning if everything is missing from the reference.
  ref_indet <- ploidies_ref == 0
  if(any(ref_indet == TRUE)) {
    message(
      "All allelic data missing for the following loci in the reference. Ploidy inferred from mixture. ",
      paste(names(ref_indet)[ref_indet], collapse = ", ")
    )
  }

  # also give a warning if everything is missing from the mixture
  mix_indet <- ploidies_mix == 0
  if(any(mix_indet == TRUE)) {
    message(
      "All allelic data missing for the following loci in the mixture. Ploidy inferred from reference. ",
      paste(names(mix_indet)[mix_indet], collapse = ", ")
    )
  }

  # now get the ploidies as 1 or 2, even if it is a 0 in either ref or mix
  ploidies <- ploidies_ref
  ploidies[ploidies == 0] <- ploidies_mix[ploidies == 0]


  ploidy_mismatch <- (ploidies_ref != ploidies_mix)[(ploidies_mix != 0) & (ploidies_ref != 0)]
  if (any(ploidy_mismatch)) {
    stop("Ploidy mismatch in reference and mixture data sets at loci ", which(ploidy_mismatch))
  }

  # check that known_collections are OK if they exist
  has_kc <- check_known_collections(reference, mixture)
  if (method == "PB" & has_kc) {
    stop("Sorry! Method PB not currently available with mixture individuals having known_collection.")
  }

  # save an untouched version of reference and gen_start_col (to be used for self-assignment)
  orig_reference <- reference
  orig_gen_start_col <- gen_start_col

  # once we are sure that repunit and collection are characters, turn them to factors
  reference$repunit <- factor(reference$repunit, levels = unique(reference$repunit))
  reference$collection <- factor(reference$collection, levels = unique(reference$collection))


  # Eric has simplified the interface.  We never expect the user to ask for
  # a trace of the individual PofZ values.  But we will always return an unthinned
  # trace of pi

  if (sample_int_Pi == 0) {
    stop("Sorry, you can't have sample_int_Pi == 0")
  }
  if ((reps + burn_in) / sample_int_Pi < 10) {
    stop("Sorry, sample_int_Pi can't be so large that you take fewer than 10 samples during burn-in and reps")
  }

  sample_int_PofZ = 0
  sample_int_omega = 0
  sample_int_rho = 0
  sample_int_PofR = 0


  # check that reference and mixture data sets have identical column names
  if (any(names(reference) != names(mixture))) stop("reference and mixture data frames differ in structure; check # of columns and variable names")

  # check to make sure that the type of each of the locus colums is the same
  type_cols_differ <- sapply(reference[-(1:(gen_start_col - 1))], class) != sapply(mixture[-(1:(gen_start_col - 1))], class)
  if (any(type_cols_differ)) stop("Data types of locus columns differ between reference and mixture at: ",
                                  paste(names(type_cols_differ[type_cols_differ]), collapse = ", "), ". Please fix that and rerun.")

  # check for a valid sampling method
  if (method != "MCMC" && method != "PB" && method != "BR") stop("invalid selection of mixture proportion estimation algorithm: please choose 'PB', 'MCMC', 'BR'")

  # get the number of missing and non-missing loci for the mixture fish and hold it
  # till the end, when we join it on there
  mix_num_loci <- count_missing_data(mixture, orig_gen_start_col)

  # Check for valid prior on pi
  if(!is.null(pi_prior)) {
    if(any(is.na(pi_prior))) stop("Custom pi_prior vectors may not contain NA values")

    if(is.numeric(pi_prior)) {
      if(length(pi_prior) != length(unique(reference$collection))) stop("Length of numeric pi prior vector does not match number of collections in reference dataset")

    } else if(is.data.frame(pi_prior)) {
      if (!("collection") %in% names(pi_prior)) stop("Missing column \"collection\" in pi_prior")
      if (!("pi_param") %in% names(pi_prior)) stop("Missing column \"pi_param\" in pi_prior")

      valid_colls <- c(as.character(unique(reference$collection)), "DEFAULT_PI")
      if(!(setequal(pi_prior$collection, intersect(pi_prior$collection, valid_colls)))) stop("Invalid collection name \"", setdiff(pi_prior$collection, intersect(pi_prior$collection, valid_colls)), "\" in pi_prior")

    } else stop("pi_prior must be NULL, a numeric vector, or a data frame")
  }


  # check to make sure pi_init is properly taken care of if non-null
  if(!is.null(pi_init)) {
    if (!("collection") %in% names(pi_init)) stop("Missing column \"collection\" in pi_init")
    if (!("pi_init") %in% names(pi_init)) stop("Missing column \"pi_init\" in pi_init")
    stopifnot(!duplicated(pi_init$collection)) # make sure no values are duplicated
    stopifnot(all(pi_init$pi_init > 0))
    all_colls <- as.character(unique(reference$collection))
    not_in_pi_init <- setdiff(pi_init$collection, all_colls)
    if(length(not_in_pi_init) > 0) stop("Missing these collections from pi_init collection column: ",
                                        paste(not_in_pi_init, collapse = ", "))
    in_pi_init_wrongly <- setdiff(all_colls, pi_init$collection)
    if(length(in_pi_init_wrongly) > 0) stop("Missing these collections from pi_init collection column: ",
                                            paste(in_pi_init_wrongly, collapse = ", "))

    named_pi_init <- pi_init$pi_init
    names(named_pi_init) <- pi_init$collection
    named_pi_init <- named_pi_init / sum(named_pi_init)
  } else {
    named_pi_init <- NULL
  }

  ## cleaning and summarizing data ##
  message("Collating data; compiling reference allele frequencies, etc.", appendLF = FALSE)

  time1 <- system.time({
    # any existing sample_type columns are removed, to be rewritten based on data frame
    if (any(names(reference) == "sample_type") || any(names(mixture) == "sample_type")) {
      reference <- dplyr::select(reference, -sample_type)
      mixture <- dplyr::select(mixture, -sample_type)
      gen_start_col <- gen_start_col - 1
    }


    # create single data frame for further processing
    D <- rbind(reference, mixture)
    sample_type <- c(rep("reference", nrow(reference)), rep("mixture", nrow(mixture)) )
    D <- cbind(sample_type, D)
    gen_start_col <- gen_start_col + 1



    # do all the cleaning and prepping necessary for inserting the reference
    # fish into the params, and grabbing a few more necessary variables
    clean <- tcf2long(D, gen_start_col)
    rac <- reference_allele_counts(clean$long)
    ac <- a_freq_list(rac)
    coll_N <- rep(0, ncol(ac[[1]])) # the number of individuals in each population; not applicable for mixture samples

    colls_by_RU <- dplyr::filter(clean$clean_short, sample_type == "reference") %>%
      droplevels() %>%
      dplyr::count(repunit, collection) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n > 0) %>% # this is now necessary for the new dplyr that counts zeroes in factor combinations
      dplyr::select(-n)


    # here we want to get a tibble of the collection names in the order in which
    # they occur in the reference once it is squashed down.  This is the levels of reference$collection
    # at this point.  And then we have to add the reporting units back on there.
    COLLS_AND_REPS_TIBBLE_CHAR <- tibble::tibble(collection = levels(reference$collection)) %>%
      dplyr::left_join(colls_by_RU %>% dplyr::mutate(repunit = as.character(repunit), collection = as.character(collection) ), by = "collection") %>%
      dplyr::select(repunit, collection)


    # COLLS_AND_REPS_TIBBLE_CHAR <- colls_by_RU %>%
    #   dplyr::mutate(repunit = as.character(repunit),
    #                 collection = as.character(collection))


    PC <- rep(0, length(unique(colls_by_RU$repunit)))
    for (i in 1:nrow(colls_by_RU)) {
      PC[colls_by_RU$repunit[i]] <- PC[colls_by_RU$repunit[i]] + 1
    }
    RU_starts <- c(0, cumsum(PC))
    RU_vec <- as.integer(colls_by_RU$collection)
    names(RU_vec) <- as.character(colls_by_RU$collection)

    # Finally, we want to prepare the alleles carried by the reference samples
    # so that we can compute the z-score statistics from them.
    ref_I <- allelic_list(cs = clean$clean_short, ac = ac, samp_type = "reference")$int

    # and get a vector of the collections those reference individuals are from
    ref_PO <- clean$clean_short %>%
      dplyr::filter(sample_type == "reference") %>%
      .$collection %>%
      as.integer()

    # and then make a params structure for doing the self-assignment to get the z-scores
    sa_params <- list_diploid_params(ac, ref_I, ref_PO, coll_N, RU_vec, RU_starts, alle_freq_prior = alle_freq_prior)
    sa_params$locus_names <- names(ac)
    sa_params$ploidies <- as.integer(unname(ploidies[sa_params$locus_names]))

  }) # close time 1 block
  message("   time: ", sprintf("%.2f", time1["elapsed"]), " seconds")



  # Now we have to get the locus-specific means and variances for computing
  # the z-scores of each individual.
  message("Computing reference locus specific means and variances for computing mixture z-scores", appendLF = FALSE)
  time_sa <- system.time({
    # ultimately, we should be able to toss the above stuff, as it has been superseded by the
    # locus-specific approach, but we will keep it around for comparison, for now...
    locus_means_and_vars <- per_locus_means_and_vars(sa_params)

  })
  message("   time: ", sprintf("%.2f", time_sa["elapsed"]), " seconds")


  #### Starting to do things to the mixture samples, now ####

  # now, we are going to break the mixture samples up into a list of data frames, each
  # named by the collection of the mixture sample, and then we are going to spew all of
  # those through the following code, and then bind them all together in the end.
  mixture_colls_list <- clean$clean_short %>%
    dplyr::filter(sample_type == "mixture") %>%
    droplevels() %>%
    split(., .$collection)



  # now, if total_catch_tib is not NULL, we break it into a named list of integer vectors
  if(!is.null(total_catch_tib)) {

    # turn total_catch_tib into a list named by collection which has recycled each
    # value up to the number of MCMC reps we are doing, and taken the ceiling of those values
    # to ensure they are integers
    tmp_list <- total_catch_tib$tot_catch
    names(tmp_list) <- total_catch_tib$collection
    tot_catch_list <- lapply(tmp_list, function(x) ceiling(rep(x, length.out = reps)))

    # then get the elements of that list in the same order as in mixture_colls_list
    tot_catch_list <- tot_catch_list[names(mixture_colls_list)]
  } else { # If total_catch_tib is NULL, then we make a tot_catch_list that has -1 for each.  This let's us test what to do once it is passed into Rcpp
    tot_catch_list <- lapply(mixture_colls_list, function(x) -1L)
  }



  # cycle over the different mixture collections and deal with each, in turn...
  big_output_list <- lapply(names(mixture_colls_list), function(name) {

    little_mix <- mixture_colls_list[[name]]
    tot_catch_vec <- tot_catch_list[[name]]

    # make some variables to handle variable_prob_is_catch == TRUE
    if(variable_prob_is_catch == TRUE) {
      doin_variable_prob_is_catch <- 1L
      prob_is_catch_vec <- little_mix$prob_is_catch
    } else {
      doin_variable_prob_is_catch <- 0L
      prob_is_catch_vec <- 1.0  # this is just a dummy to put in there so the Rcpp function does not choke
    }

    message("Working on mixture collection: ", little_mix$collection[1], " with ", nrow(little_mix), " individuals")
    if(!is.null(total_catch_tib)) {
      message("tot_catch_vec: ", paste(tot_catch_vec[1:12], collapse = " "), "...")
    }

    mix_I <- allelic_list(little_mix, ac, samp_type = "mixture")$int
    coll <- rep(0,length(mix_I[[1]]$a))  # populations of each individual in mix_I; not applicable for mixture samples



    # while we are at it, store the names of the Mixture individuals and the collections and repunits
    MIXTURE_INDIV_TIBBLE <- tibble::tibble(indiv = as.character(little_mix$indiv))



    params <- list_diploid_params(ac, mix_I, coll, coll_N, RU_vec, RU_starts, alle_freq_prior = alle_freq_prior)
    params$locus_names <- names(ac)
    params$ploidies <- as.integer(unname(ploidies[params$locus_names]))

    ## create priors on pi, if no non-default priors are submitted
    if(is.null(pi_prior)) {
      lambda <- rep(pi_prior_sum / params$C, params$C)
    } else {
      if(is.numeric(pi_prior)) {
        lambda <- pi_prior
      } else {
        lambda <- custom_pi_prior(P = pi_prior, C = data.frame(collection = colnames(ac[[1]])))
      }
    }

    ## calculate genotype log-Likelihoods for the mixture individuals ##
    message("  calculating log-likelihoods of the mixture individuals.", appendLF = FALSE)
    time2 <- system.time({
      logl <- geno_logL(params)

      # we have to be a little careful about making the scaled likelihoods, because we can
      # run into some underflow issues.  Actually that isn't really true, I think.  What
      # I was doing before was actually causing overflow.  Now I will just sweep the maxes
      # out, rather than the means.
      logl_col_maxes <- apply(logl, 2, max)
      logl_swept <- sweep(logl, 2, logl_col_maxes)

      SL <- apply(exp(logl_swept), 2, function(x) x/sum(x))

      # here we sum up the means and variances of the genotype likelihoods in each collection
      # for the loci that are non-missing in each individual, and return a C x N matrix for the meansum
      # and a C x N matrix for the varsum, in a list
      mv_sums <- rcpp_indiv_specific_logl_means_and_vars(params, locus_means_and_vars)

      # We are going to want to attach the raw likelihoods to the output for each fish
      # in the mixture.  (For computing z-scores, etc.).  So, we will make a tibble here
      # that has all of that, so we can join it back on there. logls is a matrix so we
      # just attach the indiv and collection on it and turn it all into a tibble
      logl_tibble <- tibble::tibble(indiv = base::rep(MIXTURE_INDIV_TIBBLE$indiv, each = nrow(logl)),
                                    collection = base::rep(COLLS_AND_REPS_TIBBLE_CHAR$collection, ncol(logl)),
                                    log_likelihood = base::as.vector(logl),
                                    expected_mean = base::as.vector(mv_sums$mean),
                                    expected_var = base::as.vector(mv_sums$var))


      # here is where we will modify the SL matrix to reflect the fish of known origin in the mixture
      # first get the levels of the collections in the reference
      CFL <- clean$clean_short %>%
        dplyr::filter(sample_type == "reference") %>%
        .$collection %>%
        droplevels() %>%
        levels()

      KC <- little_mix[["known_collection"]]

      SL <- modify_scaled_likelihoods_for_known_mixture_fish(SL = SL, KC = KC, CFL = CFL)

    })
    message("   time: ", sprintf("%.2f", time2["elapsed"]), " seconds")



    ## If the method is PB or MCMC, you are going to run the conditional MCMC once, at least ##
    if(method == "PB" || method == "MCMC" || (method == "BR" && !is.null(prelim_reps))) {

      # deal with initializing pi
      if(!is.null(named_pi_init)) {
        pi_init_to_use <- named_pi_init[colnames(ac[[1]])]
      } else {
        pi_init_to_use <- rep(1 / params$C, params$C)
      }

      if(method == "BR") {
        message("    performing ", prelim_reps, " initial sweeps, ", prelim_burn_in, " of which are burn-in and will not be used in computing averages to initialize starting point for method \"BR\".", appendLF = FALSE)

        time_mcmc1 <- system.time({
          out <- gsi_mcmc_1(SL = SL,
                            Pi_init = pi_init_to_use,
                            lambda = lambda,
                            reps = prelim_reps,
                            burn_in = prelim_burn_in,
                            sample_int_Pi = sample_int_Pi,
                            sample_int_PofZ = sample_int_PofZ)
        })
      } else {
        message("  performing ", reps, " total sweeps, ", burn_in, " of which are burn-in and will not be used in computing averages in method \"MCMC\"", appendLF = FALSE)
        time_mcmc1 <- system.time({
          out <- gsi_mcmc_1(SL = SL,
                            Pi_init = pi_init_to_use,
                            lambda = lambda,
                            reps = reps,
                            burn_in = burn_in,
                            sample_int_Pi = sample_int_Pi,
                            sample_int_PofZ = sample_int_PofZ,
                            sample_total_catch = as.integer(!is.null(total_catch_tib)),
                            total_catch_vals = tot_catch_vec,
                            variable_prob_is_catch = doin_variable_prob_is_catch,
                            prob_is_catch_vec = prob_is_catch_vec
                            )
        })
      }


      message("   time: ", sprintf("%.2f", time_mcmc1["elapsed"]), " seconds")
    }
    ## If the method is BR, you are going to run the MCMC with baseline resampling
    if (method == "BR") {

      message("  performing ", reps, " sweeps of method \"BR\", ",  burn_in, " sweeps of which are burn-in.", appendLF = FALSE)

      # deal with initializing pi
      if(!is.null(prelim_reps)) {
        pi_init_to_use <- out$mean$pi
      } else if(!is.null(named_pi_init)) {
        pi_init_to_use <- named_pi_init[colnames(ac[[1]])]
      } else {
        pi_init_to_use <- rep(1 / params$C, params$C)
      }

      time_mcmc2 <- system.time({
        out <- gsi_mcmc_fb(par_list = params,
                          Pi_init = pi_init_to_use,
                          lambda = lambda,
                          reps = reps,
                          burn_in = burn_in,
                          sample_int_Pi = sample_int_Pi,
                          sample_int_PofZ = sample_int_PofZ)
        names(out) <- c("mean","sd","trace")
      })
      message("   time: ", sprintf("%.2f", time_mcmc2["elapsed"]), " seconds")
    }



    ## block of code for estimating mixture using parametric bootstrap ##
    if (method == "PB") {

      # get the reporting unit proportion estimates from the original MCMC
      pi_mcmc <- out$mean$pi
      rho_mcmc <- lapply(1:(length(params$RU_starts) - 1), function(ru){
        sum(pi_mcmc[params$RU_vec[(params$RU_starts[ru] + 1):params$RU_starts[ru + 1]]])
      }) %>% unlist()



      message("  performing ", pb_iter, " bootstrapping rounds for method \"PB\"", appendLF = FALSE)
      time_pb <- system.time({
        # we have to pull the training " a" and " b"'s off the locus names
        # that got stuck there by tcf2long.  In fact, we have to rename the loci as they are in D
        little_mix_forpb <- little_mix
        names(little_mix_forpb)[gen_start_col:ncol(little_mix_forpb)] <- names(D)[gen_start_col:ncol(D)]

        ref_tmp <- D %>%
          dplyr::filter(sample_type == "reference") %>%
          dplyr::mutate(repunit = as.character(repunit), collection = as.character(collection))

        mix_tmp <- little_mix_forpb %>%
          dplyr::mutate(repunit = as.character(repunit), collection = as.character(collection))

        bootD <- rbind(ref_tmp, mix_tmp)  # bung reference and small mixture data together as required for bootstrap_rho
        boot_out <- bootstrap_rho(rho_est = rho_mcmc,
                                  pi_est = pi_mcmc,
                                  D = bootD,
                                  gen_start_col = gen_start_col,
                                  niter = pb_iter,
                                  reps = reps,
                                  burn_in = burn_in,
                                  pi_prior = pi_prior,
                                  pi_prior_sum = pi_prior_sum)

        out$mean$bootstrap_rho <- boot_out
      })
      message("   time: ", sprintf("%.2f", time_pb["elapsed"]), " seconds")
    }

    message("  tidying output into a tibble.", appendLF = FALSE)
    time_tidy <- system.time({
      ## Now for all methods we tidy up the out variable ##
      # get a tidy pi data frame #
      pi_tidy <- tidy_mcmc_coll_rep_stuff(field = out$mean,
                                          p = "pi",
                                          pname = "pi",
                                          car_tib = COLLS_AND_REPS_TIBBLE_CHAR)


      # then get a tidy PofZ
      pofz_tidy <- tidy_mcmc_pofz(input = out$mean$PofZ,
                                  pname = "PofZ",
                                  car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
                                  mix_indiv_tib = MIXTURE_INDIV_TIBBLE) %>%
        dplyr::left_join(logl_tibble, by = c("collection", "indiv"))

      # and a tidy trace of the Pi vectors
      traces_tidy <- tidy_pi_traces(input = out$trace$pi,
                                    pname = "pi",
                                    car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
                                    interval = sample_int_Pi)

      # and, if we did sampling of stock-specific total catch we get that here too, along with
      # the posterior-predictive remaining catch and the current allocations
      if(!is.null(total_catch_tib)) {
        sstc_traces <- tidy_pi_traces(
          input = out$trace$SSTC,
          pname = "SSTC",
          car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
          interval = 1
        )
        pprc_traces <- tidy_pi_traces(
          input = out$trace$PPRC,
          pname = "PPRC",
          car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
          interval = 1
        )
        allocation_count_traces <- tidy_pi_traces(
          input = out$trace$CA,
          pname = "CA",
          car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
          interval = 1
        )
      }

      ## if it was BR, we have further tidying to do to add the updated allele frequencies
      if (method == "BR") {
        if(sample_theta == TRUE) {
          theta_tidy <- lapply(1:params$L, FUN = function(x) {
            locfreqs <- tibble::as_tibble(matrix(out$mean$theta[(params$CA[x]*params$C+1):(params$CA[x+1]*params$C)],
                                                 nrow = params$A[x], ncol = params$C))
            names(locfreqs) <- colnames(ac[[x]])
            locfreqs$allele <- rownames(ac[[x]])
            locfreqs <- dplyr::select(locfreqs, allele, dplyr::everything())
          }) %>%
            dplyr::bind_rows()
          theta_tidy$locus <- rep(names(ac), times = params$A)
          theta_tidy <- dplyr::select(theta_tidy, locus, dplyr::everything())
        } else {
          theta_tidy <- NULL
        }
      }

      ## and if it was PB, we instead tidy up the bootstrap_rhos ##
      bootstrap_rhos <- NULL
      if (method == "PB") {
        bootstrap_rhos <- tibble::tibble(repunit = unique(COLLS_AND_REPS_TIBBLE_CHAR$repunit),
                                         bs_corrected_repunit_ppn = out$mean$bootstrap_rho)
      }


    })
    message("   time: ", sprintf("%.2f", time_tidy["elapsed"]), " seconds")

    # in the end, send back a list of these things
    if(method == "MCMC" || method == "PB") {
      ret <- list(mixing_proportions = pi_tidy,
           indiv_posteriors = pofz_tidy,
           mix_prop_traces = traces_tidy,
           bootstrapped_proportions = bootstrap_rhos)
    } else {
      ret <- list(mixing_proportions = pi_tidy,
           indiv_posteriors = pofz_tidy,
           mix_prop_traces = traces_tidy,
           allele_frequencies = theta_tidy)
    }

    # and, if we did stock-specific total catch sampling, send that back too
    # along with its component parts
    if(!is.null(total_catch_tib)) {
      ret$stock_specific_total_catch_traces <- sstc_traces
      ret$posterior_predictive_remaining_catch_traces <- pprc_traces
      ret$allocation_count_traces <- allocation_count_traces
    }
    ret   # return the list
  })
  # be sure to get some names onto that list
  names(big_output_list) <- names(mixture_colls_list)



  # phew.  At the end of that, we are going to bind_rows so everything is tidy, and we
  # add missing data numbers to the indiv_posteriors, too.  And we compute the z-score
  # for each individual too...
  ret <- list(
    mixing_proportions = lapply(big_output_list, function(x) x$mixing_proportions) %>%
      dplyr::bind_rows(.id = "mixture_collection"),

    indiv_posteriors = lapply(big_output_list, function(x) x$indiv_posteriors) %>%
      dplyr::bind_rows(.id = "mixture_collection") %>%
      dplyr::mutate(z_score = (log_likelihood - expected_mean) / sqrt(expected_var)) %>%
      dplyr::select(-expected_mean, -expected_var) %>%
      dplyr::left_join(., mix_num_loci, by = "indiv"),

    mix_prop_traces = lapply(big_output_list, function(x) x$mix_prop_traces) %>%
      dplyr::bind_rows(.id = "mixture_collection")
  )

  if(method == "MCMC" || method == "PB") {
    ret$bootstrapped_proportions <- lapply(big_output_list, function(x) x$bootstrapped_proportions) %>%
      dplyr::bind_rows(.id = "mixture_collection")
  } else {
    ret$allele_frequencies = lapply(big_output_list, function(x) x$allele_frequencies) %>%
      dplyr::bind_rows(.id = "mixture_collection")
  }


  # finally, here we bung together the stock-specific total catch
  if(!is.null(total_catch_tib)) {
    ret$stock_specific_total_catch_traces <- lapply(big_output_list, function(x) x$stock_specific_total_catch_traces) %>%
      dplyr::bind_rows(.id = "mixture_collection")
    ret$posterior_predictive_remaining_catch_traces <- lapply(big_output_list, function(x) x$posterior_predictive_remaining_catch_traces) %>%
      dplyr::bind_rows(.id = "mixture_collection")
    ret$allocation_count_traces <- lapply(big_output_list, function(x) x$allocation_count_traces) %>%
      dplyr::bind_rows(.id = "mixture_collection")
  }


  ret
}

