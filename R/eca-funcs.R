

#' A helper function to check that the input data frame is OK
#'
#' Just checks to make sure that column types are correct
#' @param D the data frame
#' @param gen_start_col  the column in which the genetic data starts
#' @param type For writing errors, supply "mixture" or "reference" as appropriate.
#' @keywords internal
check_refmix <- function(D, gen_start_col, type = "reference") {

  # first check to make sure that the repunit, collection, and indiv columns are present
  if(!("repunit") %in% names(D)) stop("Missing column \"repunit\" in", type)
  if(!("collection") %in% names(D)) stop("Missing column \"collection\" in", type)
  if(!("indiv") %in% names(D)) stop("Missing column \"indiv\" in", type)

  # now check to see if any of those are not character vectors
  if(!is.character(D$repunit)) stop("Column \"repunit\" must be a character vector.  It is not in ", type, " data frame")
  if(!is.character(D$collection)) stop("Column \"collection\" must be a character vector.  It is not in ", type, " data frame")
  if(!is.character(D$indiv)) stop("Column \"indiv\" must be a character vector.  It is not in ", type, " data frame")

  # now, check to make sure that all the locus columns are character or integer:
  tmp <- D[, -(1:(gen_start_col - 1))]
  char_or_int <- sapply(tmp, is.character) | sapply(tmp, is.integer)
  if(any(!char_or_int)) {
    stop("All locus columns must be of characters or integers.  These in ", type, " are not: ",
         paste(names(char_or_int[!char_or_int]), collapse = ", "))
  }

}


#' A helper function to tidy up the output from the gsi_mcmc functions
#'
#' This makes a tidy data frame of stuff, and also changes things back to
#' factors, if the levels are provided.
#' @param field  The output to tidy (i.e.. out$mean)
#' @param p the name of the parameter whose values you want to extract (like "pi")
#' @param pname the name that you want the parameter to be called in the output
#' @param car_tib  a tibble with repunit and collection in the order they appear in the output
#' @param coll_levs a vector of levels for the collection (or NULL if none)
#' @param repu_levs a vector of levels for the repunit (or NULL if none)
#' @keywords internal
tidy_mcmc_coll_rep_stuff <- function(field, p, pname, car_tib, coll_levs = NULL, repu_levs = NULL) {
  ret <- tibble::tibble(collection = car_tib$collection, value = field[[p]]) %>%
    dplyr::left_join(car_tib, ., by = "collection")

  # change the name
  names(ret)[names(ret) == "value"] <- pname

  if (!is.null(coll_levs)) {
    ret$collection = factor(ret$collection, levels = coll_levs)
  }
  if (!is.null(repu_levs)) {
    ret$repunit = factor(ret$repunit, levels = repu_levs)

  }
  ret
}






#' A helper function to tidy up the PofZ-like output from the gsi_mcmc functions
#'
#' This makes a tidy data frame of stuff, and also changes things back to
#' factors, if the levels are provided.
#' @param input  The output to tidy (i.e.. out$mean$PofZ)
#' @param pname the name that you want the parameter to be called in the output
#' @param car_tib  a tibble with repunit and collection in the order they appear in the output
#' @param mix_indiv_tib  a tibble with the individuals in the order they appear in the output
#' @param coll_levs a vector of levels for the collection (or NULL if none)
#' @param repu_levs a vector of levels for the repunit (or NULL if none)
#' @keywords internal
tidy_mcmc_pofz <- function(input, pname, car_tib, mix_indiv_tib, coll_levs = NULL, repu_levs = NULL) {
  pofz_mat <- t(input)
  colnames(pofz_mat) <- car_tib$collection

  ret <- dplyr::bind_cols(mix_indiv_tib,
                                tibble::as_tibble(pofz_mat)) %>%
    tidyr::gather(data = ., key = "collection", value = "pofz", -indiv) %>%
    dplyr::left_join(car_tib, by = "collection") %>%
    dplyr::select(indiv, repunit, collection, pofz) %>%
    dplyr::mutate(repunit = factor(repunit, levels = unique(car_tib$repunit)),   # this is heinous uglyness to get populations sorted in the order they came in the data set if they aren't factors to begin with
                  collection = factor(collection, levels = car_tib$collection)) %>%
    dplyr::arrange(indiv, repunit, collection) %>%
    dplyr::mutate(repunit = as.character(repunit),
                  collection = as.character(collection)) %>%
    dplyr::left_join(mix_indiv_tib, ., by = "indiv")

  names(ret)[names(ret) == "pofz"] <- pname

  # deal with the factor levels if present
  if (!is.null(coll_levs)) {
    ret$collection = factor(ret$collection, levels = coll_levs)
  }
  if (!is.null(repu_levs)) {
    ret$repunit = factor(ret$repunit, levels = repu_levs)

  }
  ret
}


#' a helper function to tidy up the pi-traces that come out of the mcmc functions
#'
#' This makes a tidy data frame of stuff, and also changes things back to
#' factors, if the levels are provided.
#' @param input  The output to tidy (i.e.. out$trace$pi)
#' @param pname the name that you want the parameter to be called in the output
#' @param car_tib  a tibble with repunit and collection in the order they appear in the output
#' @param mix_indiv_tib  a tibble with the individuals in the order they appear in the output
#' @param coll_levs a vector of levels for the collection (or NULL if none)
#' @param repu_levs a vector of levels for the repunit (or NULL if none)
#' @param interval the thinning interval that was used
#' @keywords internal
tidy_pi_traces <- function(input, pname, car_tib, coll_levs, repu_levs, interval) {
  ret <- lapply(input, function(x) tibble::tibble(collection = car_tib$collection,
                                            pi = x)) %>%
    dplyr::bind_rows(.id = "sweep") %>%
    dplyr::mutate(sweep = as.integer(sweep) * as.integer(interval)) %>%
    dplyr::left_join(., car_tib, by = "collection") %>%
    dplyr::select(sweep, repunit, collection, pi)

  names(ret)[names(ret) == "pi"] <- pname

  if (!is.null(coll_levs)) {
    ret$collection = factor(ret$collection, levels = coll_levs)
  }
  if (!is.null(repu_levs)) {
    ret$repunit = factor(ret$repunit, levels = repu_levs)
  }
  ret
}



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
  if(any(type_cols_differ)) stop("Data types of locus columns differ between reference and mixture at: ",
                                 paste(names(type_cols_differ[type_cols_differ]), collapse = ", "), ". Please fix that and rerun.")

  # check for a valid sampling method
  if (method != "MCMC" && method != "PB") stop("invalid selection of mixture proportion estimation algorithm: please choose 'PB', 'MCMC'")

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

    COLLS_AND_REPS_TIBBLE_CHAR <- colls_by_RU %>%
      dplyr::mutate(repunit = as.character(repunit),
                    collection = as.character(collection))


    PC <- rep(0, length(unique(colls_by_RU$repunit)))
    for (i in 1:nrow(colls_by_RU)) {
      PC[colls_by_RU$repunit[i]] <- PC[colls_by_RU$repunit[i]] + 1
    }
    RU_starts <- c(0, cumsum(PC))
    RU_vec <- as.integer(factor(colls_by_RU$collection,
                                levels = unique(colls_by_RU$collection)))

  }) # close time 1 block
  message("   time: ", sprintf("%.2f", time1["elapsed"]), " seconds")

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
      SL <- apply(exp(logl), 2, function(x) x/sum(x))
    })
    message("   time: ", sprintf("%.2f", time2["elapsed"]), " seconds")


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


    ## Now for both PB and MCMC we tidy up the out variable ##
    # get a tidy pi data frame #
    pi_tidy <- tidy_mcmc_coll_rep_stuff(field = out$mean,
                                        p = "pi",
                                        pname = "pi",
                                        car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
                                        coll_levs = NULL,
                                        repu_levs = NULL)


    # then get a tidy PofZ
    pofz_tidy <- tidy_mcmc_pofz(input = out$mean$PofZ,
                                pname = "PofZ",
                                car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
                                mix_indiv_tib = MIXTURE_INDIV_TIBBLE,
                                coll_levs = NULL,
                                repu_levs = NULL)

    # and a tidy trace of the Pi vectors
    traces_tidy <- tidy_pi_traces(input = out$trace$pi,
                                  pname = "pi",
                                  car_tib = COLLS_AND_REPS_TIBBLE_CHAR,
                                  coll_levs = NULL,
                                  repu_levs = NULL,
                                  interval = sample_int_Pi)

    ## and if it was PB, we have further tidying to do to add the bootstrap_rhos ##
    bootstrap_rhos <- NULL
    if (method == "PB") {
      bootstrap_rhos <- tibble::tibble(repunit = unique(COLLS_AND_REPS_TIBBLE_CHAR$repunit),
                                       bs_corrected_repunit_ppn = out$mean$bootstrap_rho)
    }


    # in the end, send back a list of these things
    list(mixing_proportions = pi_tidy,
         indiv_posteriors = pofz_tidy,
         mix_prop_traces = traces_tidy,
         bootstrapped_proportions = bootstrap_rhos)
  })


  # phew.  At the end of that, we are going to bind_rows so everything is tidy
  ret <- list(
    mixing_proportions = lapply(big_output_list, function(x) x$mixing_proportions) %>%
      dplyr::bind_rows(.id = "mixture_collection"),
    indiv_posteriors = lapply(big_output_list, function(x) x$indiv_posteriors) %>%
      dplyr::bind_rows(.id = "mixture_collection"),
    mix_prop_traces = lapply(big_output_list, function(x) x$mix_prop_traces) %>%
      dplyr::bind_rows(.id = "mixture_collection"),
    bootstrapped_proportions = lapply(big_output_list, function(x) x$bootstrapped_proportions) %>%
      dplyr::bind_rows(.id = "mixture_collection")
  )

  ret
}


#' Generate a samples for a mixture.
#'
#' Creates random reporting unit (rho) and collection (omega) proportions, and a
#' \code{sim_colls} vector for simulation of individual genotypes, based on the methods
#' used in Hasselman \emph{et al.} (2015)
#'
#' @param RU_starts a vector delineating the reporting units in \code{RU_vec};
#' generated by \code{tcf2param_list}
#' @param RU_vec a vector of collection indices, grouped by reporting unit;
#' generated by \code{tcf2param_list}
#' @param size The number of individuals desired in the mixture sample.  Default = 100.
#' @param alpha_repunit The dirichlet parameter for simulating the proportions of reporting units. Default = 1.5
#' @param alpha_collection The dirichlet parameter for simulating proportions of collections within reporting units. Default = 1.5
#'
#' @return a list with three elements.
#' The first two are a rho vector and an omega vector, respectively. The third is a vector of origins for
#' simulated individuals, sampled from the collections with probabilities = omega
#' @export
simulate_random_samples <- function(RU_starts, RU_vec, size = 100, alpha_repunit = 1.5, alpha_collection = 1.5) {
  rho <- gtools::rdirichlet(1, rep(alpha_repunit, length(RU_starts) - 1))
  omega <- numeric(length(RU_vec))
  omega <- sapply(1:(length(RU_starts)-1), function(x) {
    omega[RU_vec[(RU_starts[x] + 1):RU_starts[x+1]]] <- gtools::rdirichlet(1, rep(alpha_collection, RU_starts[x + 1] - RU_starts[x])) * rho[x]
  }) %>%
    cbind() %>%
    unlist()
  sim_coll = sample(RU_vec, size = size, replace = T, prob = omega)

  # finally return it all in a list
  list(rho = rho, omega = omega, sim_coll = sim_coll)
}

#' Do leave-one-out self-assignment of individuals in a reference baseline
#'
#' Returns a tidy data frame
#' @inheritParams assess_reference_loo
#' @return a tibble ...
#' @export
self_assign <- function(reference, gen_start_col) {
  # get the necessary parameters from the reference data
  params <- tcf2param_list(reference, gen_start_col, summ = T)

  # get the log-likelihoods
  logl <- t(geno_logL(par_list = params))

  # put the collection names at the top of them. To do this, we put RU_vec into sorted
  # order and then grab the names off it
  colnames(logl) <- names(sort(params$RU_vec))

  # then make a tibble of it and put the meta data (indiv, collection, repuunit) from
  # "reference" back on the results, and the gather the log-likelihoods into two columns
  # named "inferred_collection" and "log_likelihood"
  result <- reference %>%
    dplyr::select(indiv, collection, repunit) %>%
    dplyr::bind_cols(., tibble::as_tibble(logl)) %>%
    tidyr::gather(data = ., key = "inferred_collection", value = "log_likelihood", -indiv, -collection, -repunit) %>%
    dplyr::arrange(indiv, desc(log_likelihood))

  # then, if collection is a factor, we make inferred_collection a factor with the same levels
  if(is.factor(result$collection)) {
    result$inferred_collection <- factor(result$inferred_collection,
                                         levels = levels(result$collection))
  }

  # and finally, we use a join to put a column on there for "inferred_repunit".
  # this ugly thing just gets a tibble that associates repunits with collections
  repu_assoc <- result %>%
    dplyr::count(collection, repunit) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup() %>%
    dplyr::rename(inferred_collection = collection,
                  inferred_repunit = repunit)

  # and this joins the inferred_repunits column on there and then
  # orders the columns in a good way, and finally adds a column of
  # scaled likelihoods for each individual, then ungroups and returns
  # the result
  result %>%
    dplyr::left_join(., repu_assoc, by = "inferred_collection") %>%
    dplyr::select(indiv:inferred_collection, inferred_repunit, log_likelihood) %>%
    dplyr::group_by(indiv) %>%
    dplyr::mutate(scaled_likelihood = exp(log_likelihood) / sum(exp(log_likelihood))) %>%
    dplyr::ungroup()
}



#' Simulate mixtures and estimate reporting group and collection proportion estimation.
#'
#' From a reference dataset, this creates a genotype-logL matrix based on
#' simulation-by-individual with randomly drawn population proportions,
#' then uses this in two different estimates of population mixture proportions:
#' maximum likelihood via EM-algorithm and posterior mean from
#' MCMC.
#'
#' This is hard-wired at the moment to do something like Hasselman et al.
#'
#' @param reference a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns, as well as a "sample_type" column that has some "reference" entries.
#' @param gen_start_col the first column of genetic data in reference
#' @param reps  number of reps to do
#' @param mixsize the number of individuals in each simulated mixture.
#' @param seed a random seed for simulations
#' @inheritParams simulate_random_samples
#' @examples
#' ale_dev <- assess_reference_loo(alewife, 17)
#'
#' @export
assess_reference_loo <- function(reference, gen_start_col, reps = 50, mixsize = 100, seed = 5,
                                          alpha_repunit = 1.5, alpha_collection = 1.5) {

  # check that reference is formatted OK
  check_refmix(reference, gen_start_col, "reference")

  # then coerce those repunit and collection to factor to prepare them for tcf2param_list
  reference$repunit <- factor(reference$repunit, levels = unique(reference$repunit))
  reference$collection <- factor(reference$collection, levels = unique(reference$collection))

  # get the necessary parameters from the reference data
  params <- tcf2param_list(reference, gen_start_col, summ = T)

  # get a data frame that has the repunits and collections
  reps_and_colls <- reference %>%
    dplyr::group_by(repunit, collection) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::select(-n)

  # set seed
  set.seed(seed)

  # generate reps simulated data sets that each include:
  #   (1) rho's ("true" simulated mixing proportions of reporting units)
  #   (2) omegas's  ("true" simulated mixing proportions of collections)
  #   (3) sim_coll's  (a vector giving the origin of each simulated individual in the mixture)
  sim_colls <- lapply(1:reps, function(x) simulate_random_samples(params$RU_starts, params$RU_vec, size = mixsize, alpha_repunit = alpha_repunit, alpha_collection = alpha_collection))

  # now extract the true values of rho and omega from that into some data frames
  true_omega_df <- lapply(sim_colls, function(x) dplyr::data_frame(collection = levels(reference$collection), omega = x$omega)) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter))
  true_rho_df <- lapply(sim_colls, function(x) dplyr::data_frame(collection = levels(reference$repunit), rho = x$rho[1,])) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter))


  # and finally, extract the true numbers of individuals from each collection into a data frame
  true_sim_nums <- lapply(sim_colls, function(x) dplyr::data_frame(collection = names(x$sim_coll))) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter)) %>%
    dplyr::group_by(iter, collection) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup()




  #### cycle over the reps data sets and get proportion estimates from each ####
  estimates <- lapply(1:reps, function(x) {
    coll_vec <- sim_colls[[x]]$sim_coll

    # sampling SLs from the reference dataset at the individual level (like Hasselman et al. 2015)
    logL <- gprob_sim_ind(params, coll_vec)  # simulate the log-likelihood matrix of all the simmed indivs
    SL <-  apply(exp(logL), 2, function(y) y/sum(y))   # turn that into scaled likelihoods


    # get the posterior mean estimates by MCMC
    pi_out <- gsi_mcmc_1(SL = SL,
                         Pi_init = rep(1 / params$C, params$C),
                         lambda = rep(1 / params$C, params$C),
                         reps = 2000,
                         burn_in = 100,
                         sample_int_Pi = 0,
                         sample_int_PofZ = 0)



    # get the MLEs by EM-algorithm
    em_out <- gsi_em_1(SL, Pi_init = rep(1 / params$C, params$C), max_iterations = 10^6,
                       tolerance = 10^-7, return_progression = FALSE)

    # put those in a data_frame
    dplyr::data_frame(collection = levels(reference$collection),
                      post_mean = pi_out$mean$pi,
                      mle = em_out$pi
                      )
  }) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter))



  #### Now, join the estimates to the truth, re-factor everything so it is in the same order, and return ####
  ret <- dplyr::left_join(true_omega_df, true_sim_nums) %>%
    dplyr::left_join(., estimates) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n),
                  collection = factor(collection, levels = levels(reps_and_colls$collection))) %>%
    dplyr::left_join(., reps_and_colls) %>%
    dplyr::select(iter, repunit, dplyr::everything())

  # coerce repunit and collection back to character
  # and return that data frame
  ret %>%
    dplyr::mutate(collection = as.character(collection),
                  repunit = as.character(repunit))
}

#' Simulate mixtures and estimate reporting group and collection proportion estimation.
#'
#' From a reference dataset, this creates a genotype-logL matrix based on
#' simulation-by-individual with randomly drawn population proportions,
#' then uses this in two different estimates of population mixture proportions:
#' maximum likelihood via EM-algorithm and posterior mean from
#' MCMC.
#'
#' This is hard-wired at the moment to do something like Hasselman et al.
#'
#' @param reference a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns, as well as a "sample_type" column that has some "reference" entries.
#' @param gen_start_col the first column of genetic data in reference
#' @param reps  number of reps to do
#' @param mixsize the number of individuals in each simulated mixture.
#' @param seed a random seed for simulations
#' @inheritParams simulate_random_samples
#' @examples
#' ale_dev <- assess_reference_mc(alewife, 17)
#'
#' @export
assess_reference_mc <- function(reference, gen_start_col, reps = 50, mixsize = 100, seed = 5,
                                 alpha_repunit = 1.5, alpha_collection = 1.5, min_remaining = 5) {

  # check that reference is formatted appropriately
  check_refmix(reference, gen_start_col, "reference")

  reference$repunit <- factor(reference$repunit, levels = unique(reference$repunit))
  reference$collection <- factor(reference$collection, levels = unique(reference$collection))

  params <- tcf2param_list(reference, gen_start_col, summ = F)

  # get a data frame that has the repunits and collections
  reps_and_colls <- reference %>%
    dplyr::select(repunit, collection) %>%
    dplyr::group_by(repunit, collection) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(coll_int = 1:length(unique(reference$collection)))

  # set seed
  set.seed(seed)

  # get the constraints on the number of individuals to be drawn during the cross-validation
  # min_remaining individuals must be left in the reference for each collection,
  # and min_remaining * (#collections) for each reporting unit
  coll_max_draw <- reps_and_colls$n - min_remaining
  ru_max_draw <- lapply(levels(reference$repunit), function(ru){
    out <- sum(coll_max_draw[reps_and_colls$coll_int[reps_and_colls$repunit == ru]])
  }) %>% unlist()

  reps_and_colls <- dplyr::select(reps_and_colls, - coll_int)

  # Get random rhos and omegas, constrained by a minimum of min_remaing
  # reference individuals per population after the draw
  # using a stick breaking model of the Dirichlet distribution
  draw_colls <- lapply(1:reps, function(x){
    rho <- numeric(length(ru_max_draw))
    omega <- numeric(length(coll_max_draw))
    rho_sum <- 0

    for(ru in 1:length(ru_max_draw)) {
        rho[ru] <- min(ru_max_draw[ru]/mixsize,
                   (1 - rho_sum) * rbeta(1, alpha_repunit, 1.5 * (length(ru_max_draw) - ru)))
        rho_sum <- rho_sum + rho[ru]
        om_sum <- 0
        c <- 1
        for(coll in (params$RU_starts[ru] + 1):params$RU_starts[ru+1]){
          omega[params$RU_vec[coll]] <- min(coll_max_draw[params$RU_vec[coll]]/mixsize,
                             (rho[ru] - om_sum) * rbeta(1, 1.5, 1.5 * (length((params$RU_starts[ru] + 1):params$RU_starts[ru+1]) - c)))
          om_sum <- om_sum + omega[params$RU_vec[coll]]
          c <- c + 1
        }
      rho[ru] <- om_sum
    }
    # The omegas should always sum to one so long as there is a reasonable reference dataset size
    # However, could sum to less than one if the proposal for the last rho/omega is rejected
    # Therefore, include the following quick guarantee:
    rho <- rho/sum(rho)
    omega <- omega/sum(omega)
    true_n <- round(omega * mixsize,0)
    names(true_n) <- levels(reference$collection)
    list(rho=rho, omega = omega, true_n = true_n)
  })

  # now extract the true values of rho and omega from that into some data frames
  true_omega_df <- lapply(draw_colls, function(x) dplyr::data_frame(collection = levels(reference$collection), omega = x$omega)) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter))
  true_rho_df <- lapply(draw_colls, function(x) dplyr::data_frame(collection = levels(reference$repunit), rho = x$rho)) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter))

  # and finally, extract the true numbers of individuals from each collection into a data frame
  true_sim_nums <- lapply(draw_colls, function(x) dplyr::data_frame(collection = levels(reference$collection), n = x$true_n)) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter))

    reps_and_colls <- reps_and_colls %>%
    dplyr::select(-n)

  #### cycle over the reps data sets, get parameters for the new reference, and get proportion estimates from each
  estimates <- lapply(1:reps, function(x) {

    # designate random indivuals as mixture samples, based on previosly chosen proportions
    mc_data <- lapply(levels(reference$collection), function(coll){
      coll_split <- reference %>%
        dplyr::filter(collection == coll)
      mix_idx <- sample(1:nrow(coll_split), draw_colls[[x]]$true_n[coll], replace = F)
      coll_split$sample_type[mix_idx] <- "mixture"
      coll_split
      }) %>%
      dplyr::bind_rows()

    #get MCMC parameters (unique to each MC draw)
    clean <- tcf2long(mc_data, gen_start_col)
    rac <- reference_allele_counts(clean$long)
    ac <- a_freq_list(rac)
    coll_N <- rep(0, ncol(ac[[1]])) # the number of individuals in each population; not applicable for mixture samples

    colls_by_RU <- dplyr::filter(clean$clean_short, sample_type == "reference") %>%
      droplevels() %>%
      dplyr::count(repunit, collection) %>%
      dplyr::select(-n) %>%
      dplyr::ungroup()

    PC <- rep(0, length(unique(colls_by_RU$repunit)))
    for (i in 1:nrow(colls_by_RU)) {
      PC[colls_by_RU$repunit[i]] <- PC[colls_by_RU$repunit[i]] + 1
    }
    RU_starts <- c(0, cumsum(PC))
    RU_vec <- as.integer(factor(colls_by_RU$collection,
                                levels = unique(colls_by_RU$collection)))

    mix_I <- allelic_list(clean$clean_short, ac, samp_type = "mixture")$int
    coll <- rep(0,length(mix_I[[1]]$a))  # populations of each individual in mix_I; not applicable for mixture samples

    mc_params <- list_diploid_params(ac, mix_I, coll, coll_N, RU_vec, RU_starts)

    logl <- geno_logL(mc_params)
    SL <- apply(exp(logl), 2, function(x) x/sum(x))

    # get the posterior mean estimates by MCMC
    pi_out <- gsi_mcmc_1(SL = SL,
                         Pi_init = rep(1 / mc_params$C, mc_params$C),
                         lambda = rep(1 / mc_params$C, mc_params$C),
                         reps = 2000,
                         burn_in = 100,
                         sample_int_Pi = 0,
                         sample_int_PofZ = 0)



    # get the MLEs by EM-algorithm
    em_out <- gsi_em_1(SL, Pi_init = rep(1 / mc_params$C, mc_params$C), max_iterations = 10^6,
                       tolerance = 10^-7, return_progression = FALSE)

    # put those in a data_frame
    dplyr::data_frame(collection = levels(reference$collection),
                      post_mean = pi_out$mean$pi,
                      mle = em_out$pi
    )

  }) %>%
    dplyr::bind_rows(.id = "iter") %>%
    dplyr::mutate(iter = as.integer(iter))

  #### Now, join the estimates to the truth and coerce factors back to characters ####
  # first off, reps_and_colls must be converted to characters
  reps_and_colls_char <- reps_and_colls %>%
    mutate(repunit = as.character(repunit),
           collection = as.character(collection))

  ret <- dplyr::left_join(true_omega_df, true_sim_nums) %>%
    dplyr::left_join(., estimates) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    dplyr::left_join(., reps_and_colls_char) %>%
    dplyr::select(iter, repunit, dplyr::everything())

  # return that data frame
  ret
}



#' Test the effects of bias corrections on a reference dataset through cross-validation
#'
#' This is a rewrite of bias_comparison().  Eric didn't want the plotting to
#' be wrapped up in a function, and wanted to return a more informative data
#' frame.
#'
#' Takes a reference two-column genetic dataset, pulls a series of random
#' "mixture" datasets with varying reporting unit proportions from this reference,
#' and compares the results of GSI through standard MCMC, misassignment-scaled MCMC,
#' and parametric-bootstrap MCMC bias correction
#'
#' The amount of bias in reporting unit proportion calculations increases with the
#' rate of missassignment between reporting units (decreases with genetic differentiation),
#' and increases as the number of collections within reporting units becomes more uneven.
#'
#' Output from the standard Bayesian MCMC method demonstrates the level of bias to be
#' expected for the input dataset; parametric bootstrapping is an empirical method
#' for the removal of any existing bias, while the misassignment-scaled MCMC is a semi-empirical
#' method based on the rate of misassignment, which takes into account both genetic differentiation
#' and uneven collection distribution.
#'
#' @param reference a two-column format genetic dataset, with a "repunit" column
#' specifying each individual's reporting unit of origin, a "collection" column
#' specifying the collection (population or time of sampling) and "indiv" providing
#' a unique name
#' @param gen_start_col the first column containing genetic data in \code{reference}.
#' All columns should be genetic format following this column, and gene copies from the
#' same locus should be adjacent
#' @param seed the random seed for simulations
#' @param nreps The number of reps to do.
#' @param mixsize The size of each simulated mixture sample.
#'
#' @return \code{bias_comparison} returns a list; the first element is
#' a list of the relevant rho values generated on each iteration of the random "mixture"
#' creation. This includes the true rho value, the standard result \code{rho_mcmc},
#' and the parametric bootstrapped \code{rho_pb}.
#'
#' The second element is a dataframe listing summary statistics for each
#' reporting unit and estimation method. \code{mse}, the mean squared error, summarizes
#' the deviation of the rho estimates from their true value, including both bias and other variance.
#' \code{mean_prop_bias} is the average ratio of residual to true value, which gives greater
#' weight to deviations at smaller values. \code{mean_bias} is simply the average residual;
#' unlike \code{mse}, this demonstrates the direction of the bias.
#'
#' @examples
#' ale_bias <- assess_bp_bias_correction(alewife, 17)
#'
#' @export
assess_bp_bias_correction <- function(reference, gen_start_col, seed = 5, nreps = 50, mixsize = 100) {

  # check that reference is formatted appropriately
  check_refmix(reference, gen_start_col, "reference")

  reference$collection <- factor(reference$collection, levels = unique(reference$collection))
  reference$repunit <- factor(reference$repunit, levels = unique(reference$repunit))
  #get a dataframe which connects each collection to its reporting unit
  repidxs <- reference %>%
    dplyr::mutate(coll_int = as.integer(factor(reference$collection, levels = unique(reference$collection)))) %>%
    dplyr::select(repunit, coll_int) %>%
    dplyr::group_by(repunit, coll_int) %>%
    dplyr::tally()

  if(is.null(reference$sample_type)){
    sample_type <- rep("reference", nrow(reference))
    reference <- cbind(sample_type, reference)
    gen_start_col <- gen_start_col + 1
  }
  # switching any NAs in repunit and collection to prevent errors
  if(any(is.na(reference$repunit))) stop("repunit values may not be NAs" )
  if(any(is.na(reference$collection))) stop("collection values may not be NAs")
  ref_params <- tcf2param_list(reference, gen_start_col, summ = F)
  set.seed(seed)

  #fifty iterations of a system for comparing reporting unit proportion methods
  rho50 <- lapply(1:nreps, function(rr) {
    message("Starting bias_comparison rep ", rr, "   ", Sys.time())
    rho <- as.vector(gtools::rdirichlet(1, rep(1.5, length(unique(reference$repunit)))))
    #split the dataset into "reference" and "mixture", with mixture having the above rho
    drawn <- mixture_draw(reference, rhos = rho, N = mixsize, min_remaining = .0005)

    drawn_repidxs <- drawn$mixture %>%
      dplyr::group_by(repunit) %>%
      dplyr::tally()
    true_n <- drawn_repidxs$n

    # get estimates of rho from standard mcmc
    pi_mcmc <- ref_and_mix_pipeline(drawn$reference, drawn$mixture, gen_start_col, method = "MCMC")$mean$pi
    rho_mcmc <- lapply(levels(reference$repunit), function(ru){
      out <- sum(pi_mcmc[repidxs$coll_int[repidxs$repunit == ru]])
    }) %>% unlist()

    message("    Done with direct estimates. Starting bootstrap-corrected estimate...", "   ", Sys.time())
    # finally, get a bootstrap-corrected rho estimate
    delin <- rbind(drawn$reference, drawn$mixture)
    rho_pb <- bootstrap_rho(rho_mcmc, pi_mcmc, delin, gen_start_col)
    message("    Done with bootstrap-corrected estimate...", "   ", Sys.time())

    out <- list("true_rho" = rho, "true_n" = true_n, "rho_mcmc" = rho_mcmc, "rho_pb" = rho_pb)


  })

  #format data, calculate summary statistics, and generate plots
  names(rho50) <- 1:nreps
  rho50x <- lapply(rho50, tibble::as_tibble) %>%
    dplyr::bind_rows(.id = "iter")

  ret <- rho50x %>%
    dplyr::mutate(repunit = rep(unique(reference$repunit), nreps)) %>%
    dplyr::mutate(iter = as.integer(iter)) %>%
    select(iter, repunit, dplyr::everything())

  return(ret)

  if(FALSE) {  # Just removing this block.  Will wrap it up in another few functions later
    rho_data <- ret %>%
      tidyr::gather(key = "method", value = "rho_est", rho_mcmc:rho_pb)

    rho_dev <- rho_data %>%
      dplyr::mutate(dev = (true_rho - rho_est)^2) %>%
      dplyr::mutate(prop_bias = (rho_est-true_rho) / true_rho) %>%
      dplyr::mutate(bias = rho_est-true_rho) %>%
      dplyr::group_by(repunit, method) %>%
      dplyr::summarise(mse = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))


    g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = rho_est, colour = repunit)) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(repunit ~ method) +
      ggplot2::geom_abline(intercept = 0, slope = 1)
    print(g)

    # the standard mean error of each method and reporting unit
    d <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = method, y = mse, fill = repunit)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~repunit)
    print(d)

    # proportional bias makes error of equal size have larger effects
    # on small values than large values
    prop.bias <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = method, y = mean_prop_bias, fill = repunit)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~repunit)
    print(prop.bias)

    # the unscaled mean bias
    m.bias <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = method, y = mean_bias, fill = repunit)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~repunit)
    print(m.bias)

    list("rho_iterations" = rho50, "rho_dev" = rho_dev)
  }
}

