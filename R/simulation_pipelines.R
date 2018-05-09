

#' Estimate mixing proportions from reference and mixture datasets
#'
#' Takes a mixture and reference dataframe of two-column genetic data, and a
#' desired method of estimation for the population mixture proportions (MCMC, PB, or BH MCMC)
#' Returns the output of the chosen estimation method
#'
#' "MCMC" estimates mixing proportions and individual posterior
#' probabilities of assignment through Markov-chain Monte Carlo,
#' while "PB" does the same with a parametric bootstrapping correction,
#' and "BH" uses the misassignment-scaled, hierarchical MCMC.
#' All methods use a uniform 1/(# collections or RUs) prior for pi/omega and rho.
#'
#' @keywords internal
#'
#' @param reference a dataframe of two-column genetic format data, proceeded by "repunit", "collection",
#' and "indiv" columns. Does not need "sample_type" column, and will be overwritten if provided
#' @param mixture a dataframe of two-column genetic format data. Must have the same structure as
#' \code{reference} dataframe, but "collection" and "repunit" columns are ignored.
#' Does not need "sample_type" column, and will be overwritten if provided
#' @param gen_start_col the first column of genetic data in both data frames
#' @param method this must be "MCMC".   "PB" and "BH" are no longer supported in this function.
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
#' reference <- small_chinook_ref
#' mixture <- small_chinook_mix
#' gen_start_col <- 5
#'
#' # this function expects things as factors.  This function is old and needs
#' # to be replaced and deprecated.
#'
#' reference$repunit <- factor(reference$repunit, levels = unique(reference$repunit))
#' reference$collection <- factor(reference$collection, levels = unique(reference$collection))
#' mixture$repunit <- factor(mixture$repunit, levels = unique(mixture$repunit))
#' mixture$collection <- factor(mixture$collection, levels = unique(mixture$collection))
#'
#' mcmc <- ref_and_mix_pipeline(reference, mixture, gen_start_col, method = "MCMC")
#'
#' @export
ref_and_mix_pipeline <- function(reference, mixture, gen_start_col, method = "MCMC", reps = 2000, burn_in = 100, sample_int_Pi = 0, sample_int_PofZ = 0, sample_int_omega = 0, sample_int_rho = 0, sample_int_PofR = 0) {

  #check that reference and mixture data sets have identical column names
  if (any(names(reference) != names(mixture))) stop("reference and mixture data frames differ in structure; check # columns and variable names")
  # check for a valid sampling method
  if (method != "MCMC") stop("invalid selection of mixture proportion estimation algorithm: must be 'MCMC'.  'PB' and 'BH' no longer supported in this function.")
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

  # infer the ploidies
  ploidies <- get_ploidy_from_frame(D[gen_start_col:ncol(D)])


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
  for (i in 1:nrow(colls_by_RU)) {
    PC[colls_by_RU$repunit[i]] <- PC[colls_by_RU$repunit[i]] + 1
  }
  RU_starts <- c(0, cumsum(PC))
  RU_vec <- as.integer(colls_by_RU$collection)
  names(RU_vec) <- as.character(colls_by_RU$collection)
  params <- list_diploid_params(ac, mix_I, coll, coll_N, RU_vec, RU_starts)
  params$locus_names <- names(ac)
  params$ploidies <- as.integer(unname(ploidies[params$locus_names]))


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
  ref_self_params$ploidies <- as.integer(unname(ploidies[params$locus_names]))

  ref_logl <- geno_logL(ref_self_params)
  ref_SL <- apply(exp(ref_logl), 2, function(x) x/sum(x))
  ref_correctassign <- avg_coll2correctRU(ref_SL,
                                          ref_self_params$coll,
                                          ref_self_params$RU_starts,
                                          ref_self_params$RU_vec)



  # estimate population parameters based on the chosen algorithm
  if (method == "MCMC") {
    out <- gsi_mcmc_1(SL = SL,
                      Pi_init = rep(1 / params$C, params$C),
                      lambda = rep(1 / params$C, params$C),
                      reps = reps,
                      burn_in = burn_in,
                      sample_int_Pi = sample_int_Pi,
                      sample_int_PofZ = sample_int_PofZ)
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
#' @keywords internal
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
  sim_coll = sample(RU_vec, size = size, replace = TRUE, prob = omega)
  out <- list(rho = rho, omega = omega, sim_coll = sim_coll)
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
#' @param pi_prior The prior to be added to the collection allocations, in order
#' to generate pseudo-count Dirichlet parameters for the simulation of a new pi vector.
#' Non-default values should be a vector of length equal to the number of populations
#' in the reference dataset. Default value of NA leads to the
#' calculation of a symmetrical prior based on \code{pi_prior_pseudo_count_sum}.
#' @param pi_prior_pseudo_count_sum total weight on default symmetrical prior for pi.
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
#' @export
#' @keywords internal
bootstrap_rho <- function(rho_est, pi_est, D, gen_start_col, niter = 100, reps = 2000, burn_in = 100, pi_prior = NA, pi_prior_pseudo_count_sum = 1) {

  # do this to get the ploidies of the loci
  ploidies <- check_refmix(D, gen_start_col)

  D$collection <- factor(D$collection, levels = unique(D$collection))
  D$repunit <- factor(D$repunit, levels = unique(D$repunit))
  ref <- dplyr::filter(D, sample_type == "reference")
  mix <- dplyr::filter(D, sample_type == "mixture")
  repidxs <- ref %>%
    dplyr::mutate(coll_int = as.integer(factor(ref$collection, levels = unique(ref$collection)))) %>%
    dplyr::select(repunit, coll_int) %>%
    dplyr::group_by(repunit, coll_int) %>%
    dplyr::tally()

  ref_star_params <- tcf2param_list(D, gen_start_col, samp_type = "reference", summ = F, ploidies = ploidies)

  ## create symmetrical priors on pi, if no non-default priors are submitted
  if(identical(pi_prior, NA)) {
    lambda <- rep(pi_prior_pseudo_count_sum / ref_star_params$C, ref_star_params$C)
  } else lambda <- pi_prior

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
                        lambda = lambda,
                        reps = reps,
                        burn_in = burn_in,
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



