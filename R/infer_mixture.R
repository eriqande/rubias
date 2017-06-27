#' Estimate mixing proportions and origin probabilities from one or several mixtures
#'
#' Takes a mixture and reference dataframe of two-column genetic data, and a
#' desired method of estimation for the population mixture proportions (MCMC, PB)
#' Returns the output of the chosen estimation method
#'
#' "MCMC" estimates mixing proportions and individual posterior
#' probabilities of assignment through Markov-chain Monte Carlo,
#' while "PB" does the same with a parametric bootstrapping correction
#' All methods use a uniform 1/(# collections or RUs) prior for the mixing proportions.
#'
#' @param reference a dataframe of two-column genetic format data, proceeded by "repunit", "collection",
#' and "indiv" columns. Does not need "sample_type" column, and will be overwritten if provided
#' @param mixture a dataframe of two-column genetic format data. Must have the same structure as
#' \code{reference} dataframe, but "collection" and "repunit" columns are ignored.
#' Does not need "sample_type" column, and will be overwritten if provided
#' @param gen_start_col the first column of genetic data in both data frames
#' @param method a choice between "MCMC", "PB" methods for estimating mixture proportions
#' @param reps the number of iterations to be performed in MCMC
#' @param burn_in how many reps to discard in the beginning of MCMC when doing the mean calculation.
#' They will still be returned in the traces if desired.
#' @param pb_iter how many bootstrapped data sets to do for bootstrap correction using method PB.  Default
#' is 100.
#' @param sample_int_Pi how many iterations between storing the mixing proportions trace. Default is 1.
#' Can't be 0. Can't be so large that fewer than 10 samples are taken from the burn in and the sweeps.
#'
#' @return Tidy data frames in a list with the following components:
#' mixing_proportions: the estimated mixing proportions of the different collections.
#' indiv_posteriors: the posterior probs of fish being from each of the collections.
#' mix_prop_traces: the traces of the mixing proportions.  Useful for computing credible intervals.
#' bootstrapped_proportions: If using method "BH" this returns the bootstrap corrected
#' reporting unit proportions.
#'
#' @examples
#' mcmc <- infer_mixture(reference = chinook,
#'                       mixture = chinook_mix,
#'                       gen_start_col = 5,
#'                       method = "MCMC")
#'
#' @export
infer_mixture <- function(reference,
                          mixture,
                          gen_start_col,
                          method = "MCMC",
                          reps = 2000,
                          burn_in = 100,
                          pb_iter = 100,
                          sample_int_Pi = 1) {

  # check that reference and mixture are OK
  check_refmix(reference, gen_start_col, "reference")
  check_refmix(mixture, gen_start_col, "mixture")


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
  if (method != "MCMC" && method != "PB") stop("invalid selection of mixture proportion estimation algorithm: please choose 'PB', 'MCMC'")

  # get the number of missing and non-missing loci for the mixture fish and hold it
  # till the end, when we join it on there
  mix_num_loci <- count_missing_data(mixture, orig_gen_start_col)


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
      dplyr::select(-n) %>%
      dplyr::ungroup()

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
    sa_params <- list_diploid_params(ac, ref_I, ref_PO, coll_N, RU_vec, RU_starts)

  }) # close time 1 block
  message("   time: ", sprintf("%.2f", time1["elapsed"]), " seconds")



  # Now we have to do the self-assignment so that we have some raw material for computing
  # the z-scores of each individual.
  message("Performing reference self-assignment for computing mixture z-scores", appendLF = FALSE)
  time_sa <- system.time({
    sa_results <- self_assign(orig_reference, orig_gen_start_col, preCompiledParams = sa_params)
    z_score_params <- compute_reference_z_scores(sa_results, returnCollectionSummary = TRUE)
  })
  message("   time: ", sprintf("%.2f", time_sa["elapsed"]), " seconds")

  # now, we are going to break the mixture samples up into a list of data frames, each
  # named by the collection of the mixture sample, and then we are going to spew all of
  # those through the following code, and then bind them all together in the end.
  mixture_colls_list <- clean$clean_short %>%
    dplyr::filter(sample_type == "mixture") %>%
    droplevels() %>%
    split(., .$collection)


  # cycle over the different mixture collections and deal with each, in turn...
  big_output_list <- lapply(mixture_colls_list, function(little_mix) {

    message("Working on mixture collection: ", little_mix$collection[1], " with ", nrow(little_mix), " individuals")

    mix_I <- allelic_list(little_mix, ac, samp_type = "mixture")$int
    coll <- rep(0,length(mix_I[[1]]$a))  # populations of each individual in mix_I; not applicable for mixture samples



    # while we are at it, store the names of the Mixture individuals and the collections and repunits
    MIXTURE_INDIV_TIBBLE <- tibble::tibble(indiv = as.character(little_mix$indiv))



    params <- list_diploid_params(ac, mix_I, coll, coll_N, RU_vec, RU_starts)


    ## calculate genotype log-Likelihoods for the mixture individuals ##
    message("  calculating log-likelihoods of the mixture individuals.", appendLF = FALSE)
    time2 <- system.time({
      logl <- geno_logL(params)

      # we have to be a little careful about making the scaled likelihoods, because we can
      # run into some underflow issues.
      logl_col_means <- colMeans(logl)
      logl_swept <- sweep(logl, 2, logl_col_means)

      SL <- apply(exp(logl_swept), 2, function(x) x/sum(x))
    })
    message("   time: ", sprintf("%.2f", time2["elapsed"]), " seconds")

    # We are going to want to attach the raw likelihoods to the output for each fish
    # in the mixture.  (For computing z-scores, etc.).  So, we will make a tibble here
    # that has all of that, so we can join it back on there. logls is a matrix so we
    # just attach the indiv and collection on it and turn it all into a tibble
    logl_tibble <- tibble::tibble(collection = base::rep(COLLS_AND_REPS_TIBBLE_CHAR$collection, ncol(logl)),
                                  indiv = base::rep(MIXTURE_INDIV_TIBBLE$indiv, each = nrow(logl)),
                                  log_likelihood = base::as.vector(logl))


    ## regardless of whether the method is PB or MCMC, you are going to run the MCMC once, at least ##
    message("  performing ", burn_in, " burn-in and ", reps, " more sweeps of method \"MCMC\"", appendLF = FALSE)
    time_mcmc1 <- system.time({
      out <- gsi_mcmc_1(SL = SL,
                        Pi_init = rep(1 / params$C, params$C),
                        lambda = rep(1 / params$C, params$C),
                        reps = reps,
                        burn_in = burn_in,
                        sample_int_Pi = sample_int_Pi,
                        sample_int_PofZ = sample_int_PofZ)
    })
    message("   time: ", sprintf("%.2f", time_mcmc1["elapsed"]), " seconds")




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
          filter(sample_type == "reference") %>%
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
                                  burn_in = burn_in)

        out$mean$bootstrap_rho <- boot_out
      })
      message("   time: ", sprintf("%.2f", time_pb["elapsed"]), " seconds")
    }

    message("  tidying output into a tibble.", appendLF = FALSE)
    time_tidy <- system.time({
      ## Now for both PB and MCMC we tidy up the out variable ##
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

      ## and if it was PB, we have further tidying to do to add the bootstrap_rhos ##
      bootstrap_rhos <- NULL
      if (method == "PB") {
        bootstrap_rhos <- tibble::tibble(repunit = unique(COLLS_AND_REPS_TIBBLE_CHAR$repunit),
                                         bs_corrected_repunit_ppn = out$mean$bootstrap_rho)
      }

    })
    message("   time: ", sprintf("%.2f", time_tidy["elapsed"]), " seconds")

    # in the end, send back a list of these things
    list(mixing_proportions = pi_tidy,
         indiv_posteriors = pofz_tidy,
         mix_prop_traces = traces_tidy,
         bootstrapped_proportions = bootstrap_rhos)

  })


  # phew.  At the end of that, we are going to bind_rows so everything is tidy, and we
  # add missing data numbers to the indiv_posteriors, too.  And we compute the z-score
  # for each individual too...
  ret <- list(
    mixing_proportions = lapply(big_output_list, function(x) x$mixing_proportions) %>%
      dplyr::bind_rows(.id = "mixture_collection"),

    indiv_posteriors = lapply(big_output_list, function(x) x$indiv_posteriors) %>%
      dplyr::bind_rows(.id = "mixture_collection") %>%
      dplyr::left_join(., mix_num_loci, by = "indiv") %>%
      dplyr::left_join(., z_score_params, by = "collection") %>%
      dplyr::mutate(z_score = (log_likelihood - (n_non_miss_loci * L_bar)) / sqrt(n_non_miss_loci * var1)) %>%
      dplyr::select(-L_bar, -var1),

    mix_prop_traces = lapply(big_output_list, function(x) x$mix_prop_traces) %>%
      dplyr::bind_rows(.id = "mixture_collection"),

    bootstrapped_proportions = lapply(big_output_list, function(x) x$bootstrapped_proportions) %>%
      dplyr::bind_rows(.id = "mixture_collection")
  )

  ret
}

