

#' Estimate mixing proportions from reference and mixture datasets
#'
#' Takes a mixture and reference dataframe of two-column genetic data, and a
#' desired method of estimate for the population mixture proportions (EM, MCMC, or BH MCMC)
#' Returns the output of the chosen estimation method
#'
#' "EM" estimates mixing proportions and individual posterior
#' probabilities of assignment through a simple expectation maximization,
#' while "MCMC" does the same with Markov-chain Monte Carlo, and "BH" uses the misassignment-scaled,
#' Hierarchial MCMC.
#'
#' @param reference a dataframe of two-column genetic format data, proceeded by "repunit", "collection",
#' and "indiv" columns. Does not need "sample_type" column, and will be overwritten if provided
#' @param mixture a dataframe of two-column genetic format data. Must have the same structure as
#' \code{reference} dataframe, but "collection" and "repunit" columns are ignored.
#' Does not need "sample_type" column, and will be overwritten if provided
#' @param gen_start_col the first column of genetic data in both data frames
#' @param method a choice between "EM", "MCMC", and "BH" methods for estimating mixture proportions
#'
#' @return \code{mix_proportion_pipeline} returns the standard output of the chosen
#' mixing proportion estimation method (always a list)
#' @examples
#' reference <- alewife[,-1]
#' mixture <- alewife[,-1]
#' gen_start_col <- 14
#' bh <- mix_proportion_pipeline(reference, mixture, gen_start_col, method = "BH")
#' mcmc <- mix_proportion_pipeline(reference, mixture, gen_start_col, method = "MCMC")
#' em <- mix_proportion_pipeline(reference, mixture, gen_start_col, method = "EM")
#'
#' @export
ref_and_mix_pipeline <- function(reference, mixture, gen_start_col, method = "MCMC") {

  #check that reference and mixture data sets have identical column names
  if(any(names(reference) != names(mixture))) stop("reference and mixture data frames differ in structure; check # columns and variable names")
  # check for a valid sampling method
  if(method != "MCMC" && method != "EM" && method != "BH") stop("invalid selection of mixture proportion estimation algorithm: please choose 'EM', 'MCMC', or 'BH'")
  # any existing sample_type columns are removed, to be rewritten based on data frame
  if(any(names(reference) == "sample_type") || any(names(mixture) == "sample_type")) {
    reference <- dplyr::select(reference, -sample_type)
    mixture <- dplyr::select(mixture, -sample_type)
    gen_start_col <- gen_start_col - 1
  }

  # create single data frame for further processing
  D <- rbind(reference, mixture)
  sample_type <- c(rep("reference", nrow(reference)), rep("mixture", nrow(mixture)) )
  D <- cbind(sample_type, D)
  gen_start_col <- gen_start_col + 1

  # clean the data, gather allele count matrices and collection/reporting unit groups from reference data,
  # then prepare other parameters based on the mixture data
  clean <- tcf2long(D, gen_start_col)
  rac <- reference_allele_counts(clean$long)
  ac <- a_freq_list(rac)
  mix_I <- allelic_list(clean$clean_short, ac, samp_type = "mixture")$int
  coll <- rep(0,length(mix_I[[1]]$a))  # populations of each individual in mix_I; not applicable for mixture samples
  coll_N <- rep(0, ncol(ac[[1]])) # the number of individuals in each population; not applicable for mixture samples
  colls_by_RU <- dplyr::filter(clean$clean_short, sample_type == "reference") %>%
    droplevels() %>%
    dplyr::count(repunit, collection) %>%
    dplyr::select(-n)
  PC <- rep(0, length(unique(colls_by_RU$repunit)))
  for(i in 1:nrow(colls_by_RU)) {
    PC[colls_by_RU$repunit[i]] <- PC[colls_by_RU$repunit[i]] + 1
  }
  RU_starts <- c(0, cumsum(PC))
  RU_vec <- as.integer(colls_by_RU$collection)
  params <- list_diploid_params(ac, mix_I, coll, coll_N, RU_vec, RU_starts)


  # calculate genotype log-Likelihoods for the mixture individuals
  logl <- geno_logL(params)
  SL <- apply(exp(logl), 2, function(x) x/sum(x))

  #and for reference individuals, then condense to an average assigment list
  ref_I <- allelic_list(clean$clean_short, ac, samp_type = "reference")$int
  ref_coll <- as.integer(factor(reference$collection))
  ref_coll_N <- dplyr::count(reference, collection) %>%
    dplyr::select(n) %>%
    simplify2array() %>%
    as.vector()
  ref_self_params <- list_diploid_params(ac, ref_I, ref_coll, ref_coll_N, RU_vec, RU_starts)
  ref_logl <- geno_logL(ref_self_params)
  ref_SL <- apply(exp(ref_logl), 2, function(x) x/sum(x))
  ref_correctassign <- avg_coll2correctRU(ref_SL,
                                          ref_self_params$coll,
                                          ref_self_params$RU_starts,
                                          ref_self_params$RU_vec)



  # estimate population parameters based on the chosen algorithm
  if(method == "EM") {
    out <- gsi_em_1(SL, Pi_init = rep(1 / params$C, params$C), max_iterations = 10^6, tolerance = 10^-7, return_progression = FALSE)
    }
  if(method == "MCMC") {
    out <- gsi_mcmc_1(SL = SL,
                      Pi_init = rep(1 / params$C, params$C),
                      lambda = rep(1 / params$C, params$C),
                      reps = 2000,
                      burn_in = 100,
                      sample_int_Pi = 0,
                      sample_int_PofZ = 0)
  }
  if(method == "BH") {
    out <- gsi_mcmc_2(SL = SL,
                      Rho_init = rep(1 / (length(params$RU_starts) - 1), length(params$RU_starts) - 1),
                      Omega_init = rep(1 / params$C, params$C),
                      lambda_rho = rep(1 / (length(params$RU_starts) - 1), length(params$RU_starts) - 1),
                      lambda_omega = rep(1 / params$C, params$C),
                      reps = 2000,
                      burn_in = 100,
                      sample_int_omega = 0,
                      sample_int_rho = 0,
                      sample_int_PofZ = 0,
                      sample_int_PofR = 0,
                      RU_starts = params$RU_starts,
                      RU_vec = params$RU_vec,
                      coll2correctRU = ref_correctassign)
  }
out
}

#' Generate a random sim_colls vector, \emph{a la} Hasselman et al. 2015
#'
#' Creates a random sim_colls vector for simulation of individual genotypes, based on the methods
#' used in Hasselman et al. 2015, and the collection-repunit relationship vectors generated from a
#' reference
#'
#' This function is designed specifically to recreate the simulations in Hasselman, to check for
#' the bias that was observed therein. Rho (reporting unit proportions) is chosen with alphas of 1.5,
#' and omega (collection proportions) chosen with the same alpha, then scaled by the
#' corresponding rho
#'
#' @param RU_starts a vector delineating the reporting units in RU_vec
#' @param RU_vec a vector of collection indices, grouped by reporting unit
#'
#' @return \code{Hasselman_sim_colls} returns a list with three elements. The first two are
#' a rho vector and an omega vector, respectively, both with alpha parameters = 1.5 The third
#' is a vector of origins for simulated individuals, sampled from the collections with probabilities
#' = omega
#' @export
Hasselman_sim_colls <- function(RU_starts, RU_vec) {
  rho <- gtools::rdirichlet(1, c(1.5, 1.5, 1.5))
  omega <- numeric(length(RU_vec))
  omega <- sapply(1:(length(RU_starts)-1), function(x) {
    omega[RU_vec[(RU_starts[x] + 1):RU_starts[x+1]]] <- gtools::rdirichlet(1, rep(1.5, RU_starts[x + 1] - RU_starts[x])) * rho[x]
  }) %>%
    cbind() %>%
    unlist()
  sim_coll = sample(RU_vec, size = 1000, replace = T, prob = omega)
  out <- list(rho = rho, omega = omega, sim_coll = sim_coll)
}

#' Generate population proportion estimates from a reference and desired collection composition
#'
#' Using a reference allele count matrix and a vector specifying from which collection to draw
#' each individual, creates a genotype-logL matrix based on simulation-by-individual,
#' then uses this in three different estimates of population mixture proportions: EM, MCMC,
#' and hierarchical MCMC
#'
#' @param reference a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns
#' @param gen_start_col the first column of genetic data in reference
#'
#' @export
mixture_simulation_pipeline <- function(reference, gen_start_col) {

  # get the necessary parameters from the reference data
  params <- tcf2param_list(reference, gen_start_col, summ = F)
  ref_logL <- geno_logL(params)
  ref_SL <- apply(exp(ref_logL), 2, function(x) x/sum(x))
  # generate 50 sets of individuals to simulate, based on the methods in Hasselman et al. 2015
  sim_colls <- lapply(1:50, function(x) Hasselman_sim_colls(params$RU_starts, params$RU_vec))
  sim_coll_list <- lapply(1:50, function(x) sim_colls[[x]]$sim_coll)
  rhos <- lapply(1:50, function(x) as.vector(sim_colls[[x]]$rho))
  rhos <- lapply(rhos, function(x) data.frame(repunit = levels(reference$repunit), rho = x)) %>%
    dplyr::bind_rows(.id = "iter")
  omegas <- lapply(1:50, function(x) sim_colls[[x]]$omega)
  omegas <- lapply(omegas, function(x) data.frame(collection = levels(reference$collection), omega = x)) %>%
    dplyr::bind_rows(.id = "iter")

  sim_coll_tabs <- lapply(sim_coll_list, function(x) as.data.frame(table(factor(levels(reference$collection)[x], levels = levels(reference$collection))))) %>%
    setNames(1:length(sim_coll_list)) %>% dplyr::bind_rows(.id = "iter")
  names(sim_coll_tabs) <- c("iter","collection","n")
  # get proportion estimates from these sets
  props50 <- lapply(1:50, function(x) {
    coll_vec <- sim_coll_list[[x]]

    # sampling SLs from the reference dataset at the individual level (like Hasselman et al. 2015)
    logL <- gprob_sim_ind(params, coll_vec)
    SL <-  apply(exp(logL), 2, function(y) y/sum(y))

    om_out <- gsi_mcmc_2(SL = SL,
                         Rho_init = rep(1 / (length(params$RU_starts) - 1),
                                        length(params$RU_starts) - 1),
                         Omega_init = rep(1 / params$C, params$C),
                         lambda_rho = rep(1 / (length(params$RU_starts) - 1),
                                          length(params$RU_starts) - 1),
                         lambda_omega = rep(1 / params$C, params$C),
                         reps = 2000,
                         burn_in = 100,
                         sample_int_omega = 0,
                         sample_int_rho = 0,
                         sample_int_PofZ = 0,
                         sample_int_PofR = 0,
                         RU_starts = params$RU_starts,
                         RU_vec = params$RU_vec,
                         coll2correctRU = avg_coll2correctRU(ref_SL, params$coll, params$RU_starts, params$RU_vec))

    pi_out <- gsi_mcmc_1(SL = SL,
                         Pi_init = rep(1 / params$C, params$C),
                         lambda = rep(1 / params$C, params$C),
                         reps = 2000,
                         burn_in = 100,
                         sample_int_Pi = 0,
                         sample_int_PofZ = 0)

    em_out <- gsi_em_1(SL, Pi_init = rep(1 / params$C, params$C), max_iterations = 10^6,
                       tolerance = 10^-7, return_progression = FALSE)

    out <- list(bh_rho = om_out$mean$rho,
                bh_om = om_out$mean$omega,
                mc_pi = pi_out$mean$pi,
                em_pi = em_out$pi)

    names(out$bh_rho) <- levels(reference$repunit)
    out$bh_rho <- sapply(1:(length(params$RU_starts)-1), function(i) {
      rep(out$bh_rho[i], params$RU_starts[i+1]- params$RU_starts[i])
    }) %>%
      unlist()
    names(out$bh_om) <- levels(reference$collection)
    names(out$mc_pi) <- levels(reference$collection)
    names(out$em_pi) <- levels(reference$collection)
    out <- data.frame(repunit = names(out$bh_rho), collection = names(out$bh_om), as.data.frame(out))
  })

  names(props50) <- 1:50
  props50 <- props50 %>% dplyr::bind_rows(.id = "iter")



  sim_data <- dplyr::left_join(props50, rhos) %>%
    dplyr::left_join(omegas) %>%
    dplyr::left_join(sim_coll_tabs)

  rho_data <- sim_data %>%
    dplyr::group_by(iter, repunit) %>%
    dplyr::summarise(true_rho = first(rho),
                     bh_rho = first(bh_rho),
                     mcmc_rho = sum(mc_pi),
                     em_rho = sum(em_pi))

  rho_data <- rho_data %>%
    tidyr::gather(key = "method", value = "rho_est", bh_rho:em_rho)

  rho_dev <- rho_data %>%
    dplyr::mutate(dev = (true_rho - rho_est)^2) %>%
    dplyr::mutate(prop_bias = (rho_est-true_rho) / true_rho) %>%
    dplyr::mutate(bias = rho_est-true_rho) %>%
    dplyr::group_by(repunit, method) %>%
    dplyr::summarise(mean_dev = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))

  g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = rho_est, colour = repunit)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(repunit ~ method) +
    ggplot2::geom_abline(intercept = 0, slope = 1)

  print(g)

  rho_dev


}
