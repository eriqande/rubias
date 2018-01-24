
#' Test the effects of the parametric bootstrap bias correction on a reference dataset through cross-validation
#'
#' This is a rewrite of bias_comparison().  Eric didn't want the plotting to
#' be wrapped up in a function, and wanted to return a more informative data
#' frame.
#'
#' Takes a reference two-column genetic dataset, pulls a series of random
#' "mixture" datasets with varying reporting unit proportions from this reference,
#' and compares the results of GSI through standard MCMC
#' vs. parametric-bootstrap MCMC bias correction
#'
#' The amount of bias in reporting unit proportion calculations increases with the
#' rate of misassignment between reporting units (decreases with genetic differentiation),
#' and increases as the number of collections within reporting units becomes more uneven.
#'
#' Output from the standard Bayesian MCMC method demonstrates the level of bias to be
#' expected for the input data set, and parametric bootstrapping is an empirical method
#' for the removal of any existing bias.
#'
#' @param reference a two-column format genetic dataset, with a "repunit" column
#' specifying each individual's reporting unit of origin, a "collection" column
#' specifying the collection (population or time of sampling) and "indiv" providing
#' a unique name
#' @param gen_start_col the first column containing genetic data in \code{reference}.
#' All columns should be genetic format following this column, and gene copies from the
#' same locus should be adjacent
#'
#' @param seed the random seed for simulations
#' @param nreps The number of reps to do.
#' @param mixsize The size of each simulated mixture sample.
#' @param alle_freq_prior a one-element named list specifying the prior to be used when
#' generating Dirichlet parameters for genotype likelihood calculations. Valid methods include
#' \code{"const"}, \code{"scaled_const"}, and \code{"empirical"}. See
#' \code{?list_diploid_params} for method details.
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
#' \dontrun{
#' ale_bias <- assess_pb_bias_correction(alewife, 17)
#' }
#'
#' @export
assess_pb_bias_correction <- function(reference, gen_start_col, seed = 5,
                                      nreps = 50, mixsize = 100,
                                      alle_freq_prior = list("const_scaled" = 1)) {

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

  if (is.null(reference$sample_type)) {
    sample_type <- rep("reference", nrow(reference))
    reference <- cbind(sample_type, reference)
    gen_start_col <- gen_start_col + 1
  }
  # switching any NAs in repunit and collection to prevent errors
  if (any(is.na(reference$repunit))) stop("repunit values may not be NAs" )
  if (any(is.na(reference$collection))) stop("collection values may not be NAs")

  ref_params <- tcf2param_list(reference, gen_start_col, summ = F, alle_freq_prior = alle_freq_prior)

  set.seed(seed)

  #fifty iterations of a system for comparing reporting unit proportion methods
  rho50 <- lapply(1:nreps, function(rr) {
    message("Starting bias_comparison rep ", rr, "   ", Sys.time())
    rho <- as.vector(gtools::rdirichlet(1, rep(1.5, length(unique(reference$repunit)))))
    #split the dataset into "reference" and "mixture", with mixture having the above rho
    drawn <- mixture_draw(reference, rhos = rho, N = mixsize, min_remaining = .0005)

    # get the true n out of that.  There is some rigramorale here to make sure that we have
    # explicit 0's in there (that is what the left_join is all about, since tally doesn't return
    # explicit 0's for missing factor levels).
    drawn_repidxs <- drawn$mixture %>%
      dplyr::group_by(repunit) %>%
      dplyr::tally() %>%
      dplyr::mutate(repunit = as.character(repunit)) %>%
      dplyr::left_join(tibble::tibble(repunit = levels(drawn$mixture$repunit)), ., by = "repunit") %>%
      dplyr::mutate(n = ifelse(is.na(n), 0, n))

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

    out
  })

  #format data, calculate summary statistics, and generate plots
  names(rho50) <- 1:nreps
  rho50x <- lapply(rho50, tibble::as_tibble) %>%
    dplyr::bind_rows(.id = "iter")

  ret <- rho50x %>%
    dplyr::mutate(repunit = rep(unique(reference$repunit), nreps)) %>%
    dplyr::mutate(iter = as.integer(iter)) %>%
    dplyr::select(iter, repunit, dplyr::everything())

  return(ret)

  if (FALSE) {  # Just removing this block.  Will wrap it up in another few functions later
    rho_data <- ret %>%
      tidyr::gather(key = "method", value = "rho_est", rho_mcmc:rho_pb)

    rho_dev <- rho_data %>%
      dplyr::mutate(dev = (true_rho - rho_est)^2) %>%
      dplyr::mutate(prop_bias = (rho_est - true_rho) / true_rho) %>%
      dplyr::mutate(bias = rho_est - true_rho) %>%
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

