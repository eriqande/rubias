

#' Estimate mixing proportions from reference and mixture datasets
#'
#' Takes a mixture and reference dataframe of two-column genetic data, and a
#' desired method of estimation for the population mixture proportions (MCMC, PB, or BH MCMC)
#' Returns the output of the chosen estimation method
#'
#' "MCMC" estimates mixing proportions and individual posterior
#' probabilities of assignment through Markov-chain Monte Carlo,
#' while "PB" does the same with a parametric bootstrapping correction,
#' and "BH" uses the misassignment-scaled, hierarchial MCMC.
#' All methods use a uniform 1/(# collections or RUs) prior for pi/omega and rho.
#'
#' @param reference a dataframe of two-column genetic format data, proceeded by "repunit", "collection",
#' and "indiv" columns. Does not need "sample_type" column, and will be overwritten if provided
#' @param mixture a dataframe of two-column genetic format data. Must have the same structure as
#' \code{reference} dataframe, but "collection" and "repunit" columns are ignored.
#' Does not need "sample_type" column, and will be overwritten if provided
#' @param gen_start_col the first column of genetic data in both data frames
#' @param method a choice between "MCMC", "PB" and "BH" methods for estimating mixture proportions
#' @param reps the number of iterations to be performed in MCMC
#' @param burn_in how many reps to discard in the beginning of MCMC when doing the mean calculation.
#' They will still be returned in the traces if desired.
#' @param sample_int_Pi the number of reps between samples being taken for pi traces. If 0
#' no traces are taken. Only used in methods "MCMC" and "PB".
#' @param sample_int_PofZ the number of reps between samples being taken for the posterior
#' traces of each individual's collection of origin. If 0 no trace samples are taken.
#' Used in all methods
#' @param sample_int_omega the number of reps between samples being taken for
#' collection proportion traces. If 0 no traces are taken. Only used in method "BH"
#' @param sample_int_rho the number of reps between samples being taken for
#' reporting unit proportion  traces. If 0 no traces are taken. Only used in method "BH"
#' @param sample_int_PofR the number of reps between samples being taken for the posterior
#' traces of each individual's reporting unit of origin. If 0 no trace samples are taken.
#' Only used in method "BH".
#'
#' @return \code{mix_proportion_pipeline} returns the standard output of the chosen
#' mixing proportion estimation method (always a list). For method "PB",
#' returns the standard MCMC results, as well as the bootstrap-corrected
#' collection proportions under \code{$mean$bootstrap}
#' @examples
#' reference <- alewife[,-1]
#' mixture <- alewife[,-1]
#' gen_start_col <- 14
#' bh <- ref_and_mix_pipeline(reference, mixture, gen_start_col, method = "BH")
#' mcmc <- ref_and_mix_pipeline(reference, mixture, gen_start_col, method = "MCMC")
#' pb <- ref_and_mix_pipeline(reference, mixture, gen_start_col, method = "PB")
#'
#' @export
ref_and_mix_pipeline <- function(reference, mixture, gen_start_col, method = "MCMC", reps = 2000, burn_in = 100, sample_int_Pi = 0, sample_int_PofZ = 0, sample_int_omega = 0, sample_int_rho = 0, sample_int_PofR = 0) {

  #check that reference and mixture data sets have identical column names
  if(any(names(reference) != names(mixture))) stop("reference and mixture data frames differ in structure; check # columns and variable names")
  # check for a valid sampling method
  if(method != "MCMC" && method != "PB" && method != "BH") stop("invalid selection of mixture proportion estimation algorithm: please choose 'EM', 'MCMC', or 'BH'")
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
  RU_vec <- as.integer(factor(colls_by_RU$collection,
                              levels = unique(colls_by_RU$collection)))
  params <- list_diploid_params(ac, mix_I, coll, coll_N, RU_vec, RU_starts)


  # calculate genotype log-Likelihoods for the mixture individuals
  logl <- geno_logL(params)
  SL <- apply(exp(logl), 2, function(x) x/sum(x))

  #and for reference individuals, then condense to an average assigment list
  ref_I <- allelic_list(clean$clean_short, ac, samp_type = "reference")$int
  ref_coll <- as.integer(factor(reference$collection,
                                levels = unique(reference$collection)))
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
  if(method == "PB") {
    pi_out <- gsi_mcmc_1(SL = SL,
                      Pi_init = rep(1 / params$C, params$C),
                      lambda = rep(1 / params$C, params$C),
                      reps = reps,
                      burn_in = burn_in,
                      sample_int_Pi = sample_int_Pi,
                      sample_int_PofZ = sample_int_PofZ)

    pi_mcmc <- pi_out$mean$pi
    rho_mcmc <- lapply(1:(length(params$RU_starts) - 1), function(ru){
      out <- sum(pi_mcmc[params$RU_vec[(params$RU_starts[ru] + 1):params$RU_starts[ru + 1]]])
    }) %>% unlist()

    dummy_mix <- dplyr::sample_n(reference, nrow(reference), replace = TRUE)
    dummy_mix$sample_type <- rep("mixture", nrow(reference))

    boot_out <- bootstrap_rho(rho_est = rho_mcmc,
                            pi_est = pi_mcmc,
                            D = D,
                            gen_start_col = gen_start_col)
    pi_out$mean$bootstrap_rho <- boot_out
    out <- pi_out
  }
  if(method == "MCMC") {
    out <- gsi_mcmc_1(SL = SL,
                      Pi_init = rep(1 / params$C, params$C),
                      lambda = rep(1 / params$C, params$C),
                      reps = reps,
                      burn_in = burn_in,
                      sample_int_Pi = sample_int_Pi,
                      sample_int_PofZ = sample_int_PofZ)
  }
  if(method == "BH") {
    out <- gsi_mcmc_2(SL = SL,
                      Rho_init = rep(1 / (length(params$RU_starts) - 1), length(params$RU_starts) - 1),
                      Omega_init = rep(1 / params$C, params$C),
                      lambda_rho = rep(1 / (length(params$RU_starts) - 1), length(params$RU_starts) - 1),
                      lambda_omega = rep(1 / params$C, params$C),
                      reps = reps,
                      burn_in = burn_in,
                      sample_int_omega = sample_int_omega,
                      sample_int_rho = sample_int_rho,
                      sample_int_PofZ = sample_int_PofZ,
                      sample_int_PofR = sample_int_PofR,
                      RU_starts = params$RU_starts,
                      RU_vec = params$RU_vec,
                      coll2correctRU = ref_correctassign)
  }
out
}

#' Generate a random population structure and mixture sample, as in
#' Hasselman \emph{et al.} 2015
#'
#' Creates random reporting unit (rho) and collection (omega) proportions, and a
#' \code{sim_colls} vector for simulation of individual genotypes, based on the methods
#' used in Hasselman \emph{et al.} (2015)
#'
#' This function is designed specifically to recreate the simulations in Hasselman
#' \emph{et al.} (2015), to check for the bias that was observed therein.
#' Rho (reporting unit proportions) is chosen with alphas of 1.5,
#' and omega (collection proportions) chosen with the same alpha, then scaled by the
#' corresponding rho.
#'
#' @param RU_starts a vector delineating the reporting units in \code{RU_vec};
#' generated by \code{tcf2param_list}
#' @param RU_vec a vector of collection indices, grouped by reporting unit;
#' generated by \code{tcf2param_list}
#'
#' @return \code{Hasselman_sim_colls} returns a list with three elements.
#' The first two are a rho vector and an omega vector, respectively,
#' both with alpha parameters = 1.5. The third is a vector of origins for
#' simulated individuals, sampled from the collections with probabilities = omega
#' @export
Hasselman_sim_colls <- function(RU_starts, RU_vec, size = 100) {
  rho <- gtools::rdirichlet(1, rep(1.5, length(RU_starts) - 1))
  omega <- numeric(length(RU_vec))
  omega <- sapply(1:(length(RU_starts)-1), function(x) {
    omega[RU_vec[(RU_starts[x] + 1):RU_starts[x+1]]] <- gtools::rdirichlet(1, rep(1.5, RU_starts[x + 1] - RU_starts[x])) * rho[x]
  }) %>%
    cbind() %>%
    unlist()
  sim_coll = sample(RU_vec, size = size, replace = T, prob = omega)
  out <- list(rho = rho, omega = omega, sim_coll = sim_coll)
}

#' Recreate the simulations of Hasselman \emph{et al.} (2015) with two alternative
#' GSI methods
#'
#' Using a reference dataset, creates a genotype-logL matrix based on
#' simulation-by-individual with randomly drawn population proportions,
#' then uses this in three different estimates of population mixture proportions:
#' MCMC, MCMC corrected with parametric bootstrapping, and hierarchical MCMC
#'
#' The reference data set is processed with \code{tcf2param_list},
#' and the average-correct-assignment matrix is calculated using
#' \code{avg_coll2correctRU}.
#'
#' 50 Hasselman-style simulated mixture samples are then created, run through MCMC, PB,
#' and hierarchical MCMC GSI, and plotted against the rho value used to simulate
#' the mixture.
#'
#' @param reference a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns
#' @param gen_start_col the first column of genetic data in reference
#' @param seed a random seed for simulations
#'
#' @examples
#' ale_dev <- Hasselman_simulation_pipeline(alewife, 15)
#'
#' @export
Hasselman_simulation_pipeline <- function(reference, gen_start_col, seed = 5) {

  # get the necessary parameters from the reference data
  reference$collection <- factor(reference$collection, levels = unique(reference$collection))
  reference$repunit <- factor(reference$repunit, levels = unique(reference$repunit))
  params <- tcf2param_list(reference, gen_start_col, summ = F)
  ref_logL <- geno_logL(params)
  ref_SL <- apply(exp(ref_logL), 2, function(x) x/sum(x))
  avg_correct <- avg_coll2correctRU(ref_SL, params$coll, params$RU_starts, params$RU_vec)
  set.seed(seed)

  # generate 50 sets of individuals to simulate, based on the methods in Hasselman et al. 2015
  sim_colls <- lapply(1:50, function(x) Hasselman_sim_colls(params$RU_starts, params$RU_vec))
  sim_coll_list <- lapply(1:50, function(x) sim_colls[[x]]$sim_coll)
  rhos <- lapply(1:50, function(x) as.vector(sim_colls[[x]]$rho))
  rhos <- lapply(rhos, function(x) data.frame(repunit = unique(reference$repunit), rho = x)) %>%
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
                         coll2correctRU = avg_correct)

    pi_out <- gsi_mcmc_1(SL = SL,
                         Pi_init = rep(1 / params$C, params$C),
                         lambda = rep(1 / params$C, params$C),
                         reps = 2000,
                         burn_in = 100,
                         sample_int_Pi = 0,
                         sample_int_PofZ = 0)

    pi_mcmc <- pi_out$mean$pi
    rho_mcmc <- lapply(1:(length(params$RU_starts) - 1), function(ru){
      out <- sum(pi_mcmc[params$RU_vec[(params$RU_starts[ru] + 1):params$RU_starts[ru + 1]]])
    }) %>% unlist()

    dummy_mix <- dplyr::sample_n(reference, 1000, replace = TRUE)
    dummy_mix$sample_type <- rep("mixture", 1000)

    pb_out <- bootstrap_rho(rho_est = rho_mcmc,
                            pi_est = pi_mcmc,
                            D = rbind(reference, dummy_mix),
                            gen_start_col = gen_start_col)

    print(x)

    out <- list(bh_rho = om_out$mean$rho,
                bh_om = om_out$mean$omega,
                mc_pi = pi_out$mean$pi,
                pb_rho = pb_out)

    names(out$bh_rho) <- levels(reference$repunit)
    out$bh_rho <- sapply(1:(length(params$RU_starts)-1), function(i) {
      rep(out$bh_rho[i], params$RU_starts[i+1]- params$RU_starts[i])
    }) %>%
      unlist()
    names(out$pb_rho) <- levels(reference$repunit)
    out$pb_rho <- sapply(1:(length(params$RU_starts)-1), function(i) {
      rep(out$pb_rho[i], params$RU_starts[i+1]- params$RU_starts[i])
    }) %>%
      unlist()
    names(out$bh_om) <- levels(reference$collection)
    names(out$mc_pi) <- levels(reference$collection)
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
                     pb_rho = first(pb_rho))

  rho_data <- rho_data %>%
    tidyr::gather(key = "method", value = "rho_est", bh_rho:pb_rho)
  rho_data$repunit <- factor(rho_data$repunit, levels = unique(reference$repunit))

  rho_dev <- rho_data %>%
    dplyr::mutate(dev = (true_rho - rho_est)^2) %>%
    dplyr::mutate(prop_bias = (rho_est-true_rho) / true_rho) %>%
    dplyr::mutate(bias = rho_est-true_rho) %>%
    dplyr::group_by(repunit, method) %>%
    dplyr::summarise(MSE = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))

  g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = rho_est, colour = repunit)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(repunit ~ method) +
    ggplot2::geom_abline(intercept = 0, slope = 1)
print(g)

  list(rho_data, rho_dev)


}

#' Perform a parametric bootstrapping correction on an estimated rho vector
#'
#' Takes an estimate of rho, and a two-column format genetic data frame
#' containing both reference and mixture data. Returns a new rho corrected by
#' parametric bootstrapping
#'
#' @param rho_est the rho value previously estimated from MCMC GSI from the
#' provided reference and mixture data
#' @param pi_est the pi value previously estimated from MCMC GSI from the
#' provided reference and mixture data
#' @param D a two-column genetic dataframe containing the reference and mixture
#' data from which \code{rho_est} was computed; with "repunit", "collection",
#' and "indiv" columns
#' @param gen_start_col the first column of genetic data in D. All columns after
#' \code{gen_start_col} must be genetic data
#'
#' In parametric bootstrapping, \code{niter} new mixture datasets are simulated by
#' individual from the reference with reporting unit proportions \code{rho_est},
#' and the mean of their MCMC GSI outputs is used to calculate an average bias.
#' This bias is subtracted from rho_est to give the output. The number of individuals
#' in each simulated bootstrap dataset is equal to the number of "mixture" individuals
#' in \code{D}.
#'
#' @return \code{bootstrap_rho} returns a new rho value, corrected by parametric
#' bootstrapping.
#'
#' @export
bootstrap_rho <- function(rho_est, pi_est, D, gen_start_col, niter = 100) {
  D$collection <- factor(D$collection, levels = unique(D$collection))
  D$repunit <- factor(D$repunit, levels = unique(D$repunit))
  ref <- dplyr::filter(D, sample_type == "reference")
  mix <- dplyr::filter(D, sample_type == "mixture")
  repidxs <- ref %>%
    dplyr::mutate(coll_int = as.integer(factor(ref$collection, levels = unique(ref$collection)))) %>%
    dplyr::select(repunit, coll_int) %>%
    dplyr::group_by(repunit, coll_int) %>%
    dplyr::tally()

  ref_star_params <- tcf2param_list(D, gen_start_col, samp_type = "reference", summ = F)
  rho_mean <- lapply(1:niter, function(rep) {
    sim_ns <- rmultinom(n = 1, size = nrow(mix), prob = pi_est)
    sim_colls <- lapply(1:length(sim_ns), function(coll){
      rep(coll, sim_ns[coll])
    }) %>%
      unlist()
    sim_inds <- gprob_sim_ind(ref_star_params, sim_colls)
    SL <- apply(exp(sim_inds), 2, function(x) x/sum(x))
    pi_pb <- gsi_mcmc_1(SL = SL,
                        Pi_init = rep(1 / ref_star_params$C, ref_star_params$C),
                        lambda = rep(1 / ref_star_params$C, ref_star_params$C),
                        reps = 2000,
                        burn_in = 100,
                        sample_int_Pi = 0,
                        sample_int_PofZ = 0)
    rho_pb <- lapply(levels(repidxs$repunit), function(ru){
      out <- sum(pi_pb$mean$pi[repidxs$coll_int[repidxs$repunit == ru]])
    }) %>% unlist()
  }) %>%
    simplify2array() %>%
    rowMeans()

  rho_pb <- rho_est - (rho_mean - rho_est)
  # Low populations can conceivably be assigned negative values with PB, so will rescale
  rho_pb[rho_pb < 0] <- 0
  rho_pb <- rho_pb/sum(rho_pb)
  rho_pb
}



#' Test the effects of bias corrections on a reference dataset through cross-validation
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
#' the misassignment-scaled \code{rho_bh}, and the parametric bootstrapped \code{rho_pb}.
#'
#' The second element is a dataframe listing summary statistics for each
#' reporting unit and estimation method. \code{mse}, the mean squared error, summarizes
#' the deviation of the rho estimates from their true value, including both bias and other variance.
#' \code{mean_prop_bias} is the average ratio of residual to true value, which gives greater
#' weight to deviations at smaller values. \code{mean_bias} is simply the average residual;
#' unlike \code{mse}, this demonstrates the direction of the bias.
#'
#' @examples
#' ale_bias <- bias_comparison(alewife, 15)
#'
#' @export
bias_comparison <- function(reference, gen_start_col, seed = 5, nreps = 50, mixsize = 100) {

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

  # get the constraints on the number of individuals to be drawn during the cross-validation
  # a minimum of 5 individuals must be left in the reference for each collection,
  # and 5*(#collections) for each reporting unit
  #coll_max_draw <- ref_params$coll_N - 5
  #ru_max_draw <- lapply(levels(reference$repunit), function(ru){
  #out <- sum(coll_max_draw[repidxs$coll_int[repidxs$repunit == ru]])
  #}) %>% unlist()

  #fifty iterations of a system for comparing reporting unit proportion methods
  rho50 <- lapply(1:nreps, function(rr) {
    #get a random rho, constrained by a minimum of 5 individuals per population after the draw
    # using a stick breaking model of the Dirichlet distribution
    #rho <- numeric(length(ru_max_draw))
    #omega <- numeric(length(coll_max_draw))
    #rho_sum <- 0
    #for(ru in 1:length(ru_max_draw)) {
    #  rho[ru] <- min(ru_max_draw[ru]/N,
    #             (1 - rho_sum) * rbeta(1, 1.5, 1.5 * (length(ru_max_draw) - ru)))
    #  rho_sum <- rho_sum + rho[ru]
    #  om_sum <- 0
    #  c <- 1
    #  for(coll in (ref_params$RU_starts[ru] + 1):ref_params$RU_starts[ru+1]){
    #    omega[ref_params$RU_vec[coll]] <- min(coll_max_draw[ref_params$RU_vec[coll]]/N,
    #                       (rho[ru] - om_sum) * rbeta(1, 1.5, 1.5 * (length((ref_params$RU_starts[ru] + 1):ref_params$RU_starts[ru+1]) - c)))
    #    om_sum <- om_sum + omega[ref_params$RU_vec[coll]]
    #    c <- c + 1
    #  }
    # encountered a bug where if the omega proposal is rejected for the
    # last collection in a reporting unit, the acceptance of the max/N
    # causes omega to sum to less than one; the following line should
    # fix this bug for all but the last collection to be chosen
    #  rho[ru] <- om_sum
    #}
    # quick fix in case the last omega to be chosen is rejected;
    # should find a better solution
    #rho <- rho/sum(rho)
    #omega <- omega/sum(omega)
    #if(!identical(all.equal(sum(omega),1), TRUE)) print(omega)
    message("Starting bias_comparison rep ", rr, "   ", Sys.time())
    rho <- as.vector(gtools::rdirichlet(1, rep(1.5, length(unique(reference$repunit)))))
    #split the dataset into "reference" and "mixture", with mixture having the above rho
    drawn <- mixture_draw(reference, rhos = rho, N = mixsize, min_remaining = .0005)
    # get estimates of rho from standard mcmc
    pi_mcmc <- ref_and_mix_pipeline(drawn$reference, drawn$mixture, gen_start_col, method = "MCMC")$mean$pi
    rho_mcmc <- lapply(levels(reference$repunit), function(ru){
      out <- sum(pi_mcmc[repidxs$coll_int[repidxs$repunit == ru]])
    }) %>% unlist()
    # and from finagled mcmc
    rho_bh <- ref_and_mix_pipeline(drawn$reference, drawn$mixture, gen_start_col, method = "BH")$mean$rho

    message("    Done with direct estimates. Starting bootstrap-corrected estimate...", "   ", Sys.time())
    # finally, get a bootstrap-corrected rho estimate
    delin <- rbind(drawn$reference, drawn$mixture)
    rho_pb <- bootstrap_rho(rho_mcmc, pi_mcmc, delin, gen_start_col)
    message("    Done with bootstrap-corrected estimate...", "   ", Sys.time())

    out <- list("true_rho" = rho, "rho_mcmc" = rho_mcmc, "rho_bh" = rho_bh, "rho_pb" = rho_pb)


  })

  #format data, calculate summary statistics, and generate plots
  names(rho50) <- 1:nreps
  rho50x <- lapply(rho50, tibble::as_tibble) %>%
    dplyr::bind_rows(.id = "iter")

  rho50x$repunit <- rep(unique(reference$repunit), nreps)

  rho_data <- rho50x %>%
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

