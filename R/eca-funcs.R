#' Estimate mixing proportions and origin probabilities from a mixture
#'
#' THIS IS A MINOR REWRITE OF ref_and_mix_pipeline. ERIC JUST WANTED TO CHANGE
#' THE OUTPUT A LITTLE BIT SO IT IS RETURNING DATA FRAMES. AND MORE FLEXIBILITY
#' IN SPECIFYING PRIORS.
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
#' @param pb_iter how many bootstrapped data sets to do for bootstrap correction using method PB.  Default
#' is 100.
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
infer_mixture <- function(reference,
                          mixture,
                          gen_start_col,
                          method = "MCMC",
                          reps = 2000,
                          burn_in = 100,
                          pb_iter = 100,
                          sample_int_Pi = 0,
                          sample_int_PofZ = 0,
                          sample_int_omega = 0,
                          sample_int_rho = 0,
                          sample_int_PofR = 0) {

  # check that repunit and population are factors in reference
  if(!is.factor(reference$repunit)) stop("repunit column in input reference must be a factor")
  if(!is.factor(reference$collection)) stop("collection column in input reference must be a factor")

  # check that reference and mixture data sets have identical column names
  if(any(names(reference) != names(mixture))) stop("reference and mixture data frames differ in structure; check # columns and variable names")

  # check for a valid sampling method
  if(method != "MCMC" && method != "PB" && method != "BH") stop("invalid selection of mixture proportion estimation algorithm: please choose 'EM', 'MCMC', or 'BH'")

  message("collating data; compiling allele frequencies, etc.", appendLF = FALSE)

  time1 <- system.time({
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
  }) # close time 1 block
  message("   time: ", sprintf("%.2f", time1["elapsed"]), " seconds")

  # calculate genotype log-Likelihoods for the mixture individuals
  message("calculating log-likelihoods of the mixture individuals.", appendLF = FALSE)
  time2 <- system.time({
    logl <- geno_logL(params)
    SL <- apply(exp(logl), 2, function(x) x/sum(x))
  })
  message("   time: ", sprintf("%.2f", time2["elapsed"]), " seconds")


  # estimate population parameters based on the chosen algorithm
  if(method == "PB") {
    message("performing ", burn_in, " burn-in and ", reps, " more sweeps in first round of method \"PB\"", appendLF = FALSE)
    time_mcmc1 <- system.time({
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
    })
    message("   time: ", sprintf("%.2f", time_mcmc1["elapsed"]), " seconds")

    #dummy_mix <- dplyr::sample_n(reference, nrow(reference), replace = TRUE)
    #dummy_mix$sample_type <- rep("mixture", nrow(reference))

    message("performing ", pb_iter, " bootstrapping rounds for method \"PB\"", appendLF = FALSE)
    time_pb <- system.time({
      boot_out <- bootstrap_rho(rho_est = rho_mcmc,
                                pi_est = pi_mcmc,
                                D = D,
                                gen_start_col = gen_start_col,
                                niter = pb_iter)
      pi_out$mean$bootstrap_rho <- boot_out
      out <- pi_out
    })
    message("   time: ", sprintf("%.2f", time_pb["elapsed"]), " seconds")
  }
  if(method == "MCMC") {
    message("performing ", burn_in, " burn-in and ", reps, " more sweeps of method \"MCMC\"", appendLF = FALSE)
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
  }
  if(method == "BH") {
    message("performing ", burn_in, " burn-in and ", reps, " more sweeps of method \"BH\"", appendLF = FALSE)
    time_mcmc2 <- system.time({
      out <- gsi_mcmc_bh(SL = SL,
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
                        RU_vec = params$RU_vec)
    })
    message("   time: ", sprintf("%.2f", time_mcmc2["elapsed"]), " seconds")
  }
  out
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
#' ale_dev <- simulate_and_assess_reference(alewife, 15)
#'
#' @export
simulate_and_assess_reference <- function(reference, gen_start_col, reps = 50, mixsize = 100, seed = 5,
                                          alpha_repunit = 1.5, alpha_collection = 1.5) {

  # get the necessary parameters from the reference data
  params <- tcf2param_list(reference, gen_start_col, summ = F)

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
assess_bp_bias_correction <- function(reference, gen_start_col, seed = 5, nreps = 50, mixsize = 100) {

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

