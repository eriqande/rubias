
#' Simulate mixtures and estimate reporting group and collection proportions.
#'
#' From a reference dataset, this creates a genotype-logL matrix based on
#' simulation-by-individual with randomly drawn population proportions,
#' then uses this in two different estimates of population mixture proportions:
#' maximum likelihood via EM-algorithm and posterior mean from
#' MCMC.
#'
#' This is hard-wired at the moment to do something like Hasselman et al. (2015)
#'
#' @param reference a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns, as well as a "sample_type" column that has some "reference" entries
#' @param gen_start_col the first column of genetic data in \code{reference}
#' @param reps  number of reps of mixture simulation and MCMC to do
#' @param mixsize the number of individuals in each simulated mixture
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
  ret <- dplyr::left_join(true_omega_df, true_sim_nums, by = c("iter", "collection")) %>%
    dplyr::left_join(., estimates, by = c("iter", "collection")) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n),
                  collection = factor(collection, levels = levels(reps_and_colls$collection))) %>%
    dplyr::left_join(., reps_and_colls, by = "collection") %>%
    dplyr::select(iter, repunit, dplyr::everything())

  # coerce repunit and collection back to character
  # and return that data frame after renaming the variables their final form
  ret %>%
    dplyr::mutate(collection = as.character(collection),
                  repunit = as.character(repunit)) %>%
    dplyr::rename(true_pi = omega,
                  post_mean_pi = post_mean,
                  mle_pi = mle)
}
