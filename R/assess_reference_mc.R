#' Partition a reference dataset and estimate reporting group and collection proportions
#'
#' From a reference dataset, this draws (without replacement) a simulated mixture
#' dataset with randomly drawn population proportions,
#' then uses this in two different estimates of population mixture proportions:
#' maximum likelihood via EM-algorithm and posterior mean from MCMC.
#'
#' This method is referred to as "Monte Carlo cross-validation".
#' The input parameters for \code{assess_reference_mc} are more restrictive than those of
#' \code{assess_reference_loo}. Rather than allowing a \emph{data.frame} to specify Dirichlet
#' parameters, proportions, or counts for specific reporting units and collections,
#' \code{assess_reference_mc} only allows vector input (default = 1.5) for \code{alpha_repunit}
#' and \code{alpha_collection}. These inputs specify the uniform Dirichlet parameters for
#' all reporting units and collections, respectively.
#'
#' For mixture proportion generation, the rho values are first drawn using a stick-breaking
#' model of the Dirichlet distribution, but with proportions capped by \code{min_remaining}.
#' Stick-breaking is then used to subdivide each reporting unit into collections. In
#' addition to the constraint that mixture sampling without replacement cannot deplete
#' the number of individuals in each collection below \code{min_remaining}, a similar
#' constraint is placed upon the number of individuals left in reporting units,
#' determined as \code{min_remaining} * (# collections in reporting unit).
#'
#' Note that this implies that the data are only truly Dirichlet distributed when no
#' rejections based on \code{min_remaining} occur. This is a reasonable certainty with
#' sufficient reference collection sizes relative to the desired mixture size.
#'
#' @param reference a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns, as well as a "sample_type" column that has some "reference" entries.
#' @param gen_start_col the first column of genetic data in reference
#' @param reps  number of reps to do
#' @param mixsize the number of individuals in each simulated mixture.
#' @param seed a random seed for simulations
#' @param alpha_repunit The dirichlet parameter for simulating the proportions of
#' reporting units. Default = 1.5
#' @param alpha_collection The dirichlet parameter for simulating proportions of collections
#' within reporting units. Default = 1.5
#' @param min_remaining the minimum number of individuals which should be conserved in
#' each reference collection during sampling without replacement to form the simulated mixture
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

  reps_and_colls <- dplyr::select(reps_and_colls, -coll_int)

  # Get random rhos and omegas, constrained by a minimum of min_remaining
  # reference individuals per population after the draw
  # using a stick breaking model of the Dirichlet distribution
  draw_colls <- lapply(1:reps, function(x){
    rho <- numeric(length(ru_max_draw))
    omega <- numeric(length(coll_max_draw))
    rho_sum <- 0

    for (ru in 1:length(ru_max_draw)) {
      rho[ru] <- min(ru_max_draw[ru]/mixsize,
                     (1 - rho_sum) * rbeta(1, alpha_repunit, 1.5 * (length(ru_max_draw) - ru)))
      rho_sum <- rho_sum + rho[ru]
      om_sum <- 0
      c <- 1
      for (coll in (params$RU_starts[ru] + 1):params$RU_starts[ru + 1]) {
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

