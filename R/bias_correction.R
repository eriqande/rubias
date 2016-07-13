
#' Test the effects of bias correction on a particular reference dataset
#'
#' Takes a reference two-column genetic dataset, pulls a series of random
#' "mixture" datasets with varying reporting unit proportions from this reference,
#' and compares the results of GSI through standard MCMC, misassignment-scaled MCMC,
#' and parametric-bootstrap MCMC bias correction
#'
#' @param reference a two-column format genetic dataset, with a "repunit" column
#' specifying each individual's reporting unit of origin, a "collection" column
#' specifying the collection (population or time of sampling) and "indiv" providing
#' a unique name
#' @param gen_start_col the first column containing genetic data in \code{reference}.
#' All columns should be genetic format following this column, and gene copies from the
#' same locus should be adjacent
#' @param seed the random seed for simulations
#'
#' @return \code{bias_comparison} returns a list; the first element is
#' a list of the relevant rho values generated on each iteration of the random "mixture"
#' creation. This includes the true rho value, the standard result \code{rho_mcmc},
#' the misassignment-scaled \code{rho_bh}, and the parametric bootstrapped \code{rho_pb}.
#' The second element is a dataframe listing summary statistics for each
#' reporting unit and estimation method
#'
#' @examples
#' ale_bias <- bias_comparison(alewife, 15)
#'
#' @export
bias_comparison <- function(reference, gen_start_col, seed = 5) {

  #get a dataframe which connects each collection to its reporting unit
  repidxs <- reference %>%
    dplyr::mutate(coll_int = as.integer(as.factor(collection))) %>%
    dplyr::select(repunit, coll_int) %>%
    dplyr::group_by(repunit, coll_int) %>%
    dplyr::tally()

  if(is.null(reference$sample_type)){
    sample_type <- rep("reference", nrow(reference))
    reference <- cbind(sample_type, reference)
  }
  gen_start_col <- gen_start_col + 1
  ref_params <- tcf2param_list(reference, gen_start_col, summ = F)
  set.seed(seed)

  #fifty iterations of a system for comparing reporting unit proportion methods
  rho50 <- lapply(1:50, function(rr) {
    #get a random rho
    rho <- gtools::rdirichlet(1, rep(1.5, length(unique(reference$repunit))))
    #split the dataset into "reference" and "mixture", with mixture having the above rho
    drawn <- mixture_draw(reference, rhos = rho, N = 100, min_remaining = .005)
    # get estimates of rho from standard mcmc
    pi_mcmc <- ref_and_mix_pipeline(drawn$reference, drawn$mixture, 15, method = "MCMC")$mean$pi
    rho_mcmc <- lapply(levels(reference$repunit), function(ru){
      out <- sum(pi_mcmc[repidxs$coll_int[repidxs$repunit == ru]])
    }) %>% unlist()
    # and from finagled mcmc
    rho_bh <- ref_and_mix_pipeline(drawn$reference, drawn$mixture, 15, method = "BH")$mean$rho

    # finally, get a bootstrap-corrected rho estimate
    delin <- rbind(drawn$reference, drawn$mixture)
    ref_star_params <- tcf2param_list(delin, 15, samp_type = "reference", summ = F)
    rho_mean <- lapply(1:100, function(rep) {
      sim_ns <- rmultinom(n = 1, size = length(pi_mcmc), prob = pi_mcmc)
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
      rho_pb <- lapply(levels(reference$repunit), function(ru){
        out <- sum(pi_pb$mean$pi[repidxs$coll_int[repidxs$repunit == ru]])
      }) %>% unlist()
    }) %>%
      simplify2array() %>%
      rowMeans()

    rho_pb <- rho_mcmc - (rho_mean - rho_mcmc)

    out <- list("true_rho" = rho, "rho_mcmc" = rho_mcmc, "rho_bh" = rho_bh, "rho_pb" = rho_pb)
  })

  #format data, calculate summary statistics, and generate plots
  names(rho50) <- 1:50
  rho50x <- rho50 %>% dplyr::bind_rows(.id = "iter")
  rho50x$repunit <- rep(unique(reference$repunit), 50)

  rho_data <- rho50x %>%
    tidyr::gather(key = "method", value = "rho_est", rho_mcmc:rho_pb)

  rho_dev <- rho_data %>%
    dplyr::mutate(dev = (true_rho - rho_est)^2) %>%
    dplyr::mutate(prop_bias = (rho_est-true_rho) / true_rho) %>%
    dplyr::mutate(bias = rho_est-true_rho) %>%
    dplyr::group_by(repunit, method) %>%
    dplyr::summarise(sme = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))


  g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = rho_est, colour = repunit)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(repunit ~ method) +
    ggplot2::geom_abline(intercept = 0, slope = 1)
  print(g)

  # the standard mean error of each method and reporting unit
  d <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = method, y = sme, fill = repunit)) +
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
