
#' Simulate mixtures and estimate reporting group and collection proportions.
#'
#' From a reference dataset, this creates a genotype-logL matrix based on
#' simulation-by-individual with randomly drawn population proportions,
#' then uses this in two different estimates of population mixture proportions:
#' maximum likelihood via EM-algorithm and posterior mean from
#' MCMC.
#'
#'
#' @param reference a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns, as well as a "sample_type" column that has some "reference" entries
#' @param gen_start_col the first column of genetic data in \code{reference}
#' @param reps  number of reps of mixture simulation and MCMC to do
#' @param mixsize the number of individuals in each simulated mixture
#' @param seed a random seed for simulations
#' @param printSummary if TRUE a summary of the reference samples will be printed to stdout.
#' @param resampling_unit what unit should be resampled.  Currently the choices are "individuals" (the default)
#' and "gene_copies".  Using "individuals" preserves missing data patterns available in the reference data set.
#' We also have "gene_copies_with_missing" capability, but it is not yet linked into this function.
#' @param alpha_repunit If a vector, this is the dirichlet parameter for simulating
#' the proportions of reporting units. Gets recycled to the number of reporting units. Default is 1.5.
#' Otherwise, this could be a two-column data frame.  The first column must be named "repunit" and the
#' second one must be one of "dirichlet", "ppn", or "cnt", according to whether you wish to
#' specify dirichlet parameters, or proportions, or exact counts, respectively, for each population.
#' If you want to make multiple simulations, pass in a list of data frames or of individual dirichlet parameters.
#' @param alpha_collection The dirichlet parameter for simulating proportions of collections within reporting units. Default = 1.5.
#' If this is a data frame then the first column must be "collection" and the second must be one of
#' "dirichlet", "ppn", "cnt", "sub_dirichlet", "sub_ppn".  If you want to provide multiple different
#' scenarios.  You can pass them in as a list.  If alpha_repunit or alpha_collection is a list with length
#' greater than 1, the shorter will be recycled.
#' @details We still need to implement storage of individual-specific results for the simulated individuals.
#' This we can do pretty easily, but wanted to get the proportions done and out there first.  Send eric some
#' email if you really need individual specific output and he will prioritize it.
#' @examples
#' ale_dev <- assess_reference_loo(alewife, 17)
#'
#' @export
assess_reference_loo <- function(reference, gen_start_col, reps = 50, mixsize = 100, seed = 5,
                                 alpha_repunit = 1.5, alpha_collection = 1.5, resampling_unit = "individual",
                                 printSummary = FALSE) {

  if (!(resampling_unit %in% c("gene_copies", "individual"))) stop("Choice ", resampling_unit, " unknown for resampling unit.")

  # check that reference is formatted OK
  check_refmix(reference, gen_start_col, "reference")

  # then coerce those repunit and collection to factor to prepare them for tcf2param_list
  reference$repunit <- factor(reference$repunit, levels = unique(reference$repunit))
  reference$collection <- factor(reference$collection, levels = unique(reference$collection))


  # get the necessary parameters from the reference data
  params <- tcf2param_list(reference, gen_start_col, summ = printSummary)

  # get a data frame that has the repunits and collections
  reps_and_colls <- reference %>%
    dplyr::group_by(repunit, collection) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::select(-n)

  # set seed
  set.seed(seed)

  # now, take care of the alpha_repunit and alpha_collection parameters.
  if (class(alpha_repunit)[1] == "list") {
    alpha_repu_list <- alpha_repunit
  } else {
    alpha_repu_list <- list(alpha_repunit)  # make it a list of length 1 if need be
  }
  if (class(alpha_collection)[1] == "list") {
    alpha_coll_list <- alpha_collection
  } else {
    alpha_coll_list <- list(alpha_collection)  # make it a list of length 1 if need be
  }

  # at this juncture we now have lists of scenarios: alpha_coll_list and alpha_repu_list
  # which we will cycle over after recyling them to the longest lenght


  # We set up names on the lists, if absent so we know which was which.
  # Note that people can name their scenarios for easy access later
  if (is.null(names(alpha_repu_list))) {
    names(alpha_repu_list) <- 1:length(alpha_repu_list)
  }
  if (is.null(names(alpha_coll_list))) {
    names(alpha_coll_list) <- 1:length(alpha_coll_list)
  }
  # now, recycle values in those lists as need be.  We will recycle reps as well, so you can pass
  # in a vector and it will cycle through that over different scenarios.
  max_list_length = max(length(alpha_coll_list), length(alpha_repu_list))
  alpha_repu_list <- rep(alpha_repu_list, length.out = max_list_length)
  alpha_coll_list <- rep(alpha_coll_list, length.out = max_list_length)
  rep_vec <- rep(reps, length.out = max_list_length)

  # now make a list for outputting all these results (though we will bind_rows them all together at the end...)
  output_list <- list()

  # and now a good ol' fashioned for loop:
  for (scenario in 1:max_list_length) {
    # get all the current values for this iteration:
    alpha_repunit <- alpha_repu_list[[scenario]]
    alpha_collection <- alpha_coll_list[[scenario]]
    reps <- rep_vec[scenario]
    repu_scenario <- names(alpha_repu_list)[scenario]
    coll_scenario <- names(alpha_coll_list)[scenario]

    # print a message about progress
    message("++++ Starting in on repunit_scenario ", repu_scenario, " with collection scenario ", coll_scenario, " ++++")

    # generate reps simulated data sets that each include:
    #   (1) rho's ("true" simulated mixing proportions of reporting units)
    #   (2) omegas's  ("true" simulated mixing proportions of collections)
    #   (3) sim_coll's  (a vector giving the origin of each simulated individual in the mixture)
    sim_colls <- lapply(1:reps, function(x)  {
      simulate_random_samples(params$RU_starts,
                              params$RU_vec,
                              size = mixsize,
                              alpha_repunit = alpha_repunit,
                              alpha_collection = alpha_collection)
    })



    # now extract the true values of rho and omega from that into some data frames
    true_omega_df <- lapply(sim_colls, function(x) tibble::enframe(x$omega, name = "collection", value = "omega")) %>%
      dplyr::bind_rows(.id = "iter") %>%
      dplyr::mutate(iter = as.integer(iter))
    true_rho_df <- lapply(sim_colls, function(x) tibble::enframe(x$rho, name = "repunit", value = "rho")) %>%
      dplyr::bind_rows(.id = "iter") %>%
      dplyr::mutate(iter = as.integer(iter))


    # and finally, extract the true numbers of individuals from each collection into a data frame
    true_sim_nums <- lapply(sim_colls, function(x) tibble::tibble(collection = names(x$sim_coll))) %>%
      dplyr::bind_rows(.id = "iter") %>%
      dplyr::mutate(iter = as.integer(iter)) %>%
      dplyr::group_by(iter, collection) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::ungroup()




    #### cycle over the reps data sets and get proportion estimates from each ####
    estimates <- lapply(1:reps, function(x) {

      message("Doing LOO simulations rep ", x, " of ", reps)

      coll_vec <- sim_colls[[x]]$sim_coll

      # sampling SLs from the reference dataset
      if (resampling_unit == "individual") {
        logL <- gprob_sim_ind(params, coll_vec)  # simulate the log-likelihood matrix of all the simmed indivs
      }
      if (resampling_unit == "gene_copies") {
        logL <- gprob_sim_gc(params, coll_vec)
      }

      # we have to be a little careful about making the scaled likelihoods, because we can
      # run into some underflow issues.
      logl_col_means <- colMeans(logL)
      logl_swept <- sweep(logL, 2, logl_col_means)

      SL <-  apply(exp(logl_swept), 2, function(y) y/sum(y))   # turn that into scaled likelihoods


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
    ret <- dplyr::left_join(true_omega_df, true_sim_nums, by = c("iter", "collection")) %>%
      dplyr::left_join(., estimates, by = c("iter", "collection")) %>%
      dplyr::mutate(n = ifelse(is.na(n), 0, n),
                    collection = factor(collection, levels = levels(reps_and_colls$collection))) %>%
      dplyr::left_join(., reps_and_colls, by = "collection") %>%
      dplyr::select(iter, repunit, dplyr::everything())

    # coerce repunit and collection back to character
    # and return that data frame after renaming the variables to their final form
    output_list[[scenario]] <- ret %>%
      dplyr::mutate(collection = as.character(collection),
                    repunit = as.character(repunit)) %>%
      dplyr::rename(true_pi = omega,
                    post_mean_pi = post_mean,
                    mle_pi = mle) %>%
      dplyr::mutate(repunit_scenario = repu_scenario,
                    collection_scenario = coll_scenario)
  }

  # in the end we bind those all together and put the scenario columns up front
  dplyr::bind_rows(output_list) %>%
    dplyr::select(repunit_scenario, collection_scenario, dplyr::everything())
}
