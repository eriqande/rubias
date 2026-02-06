#' Make a plot comparing PofZs and Pi values for all collections (or repunits)
#'
#' @param X The return object from infer_mixture
#' @param mix_coll  character name of the mixture_collection within X
#' you want to plot results for.
#' @param by_repunit logical.  If true, do by reporting units rather than collections.
#' @param burn_in  discards the first burn_in sweeps when handling the pi_traces
#' (for the boxplots of the pi)
#' @param plot_flavor a string giving the possibilities: "normal": the there is
#' no scaling on the bar plots.  "bars_fully_expanded": the horizontal barplots (and boxplots)
#' are drawn on an expanded scale so that the larger of the sum of posteriors or the
#' sum of scaled likelihoods is set to the value of N * pi_big, where pi_big is the
#' largest mixing proportion and N is the sample size.  This allows us to see things
#' that are otherwise more squashed down for the rare stocks. "bars_expanded_and_ticked"
#' bars are expanded, but any group with a sum over individuals of PofZ or scaled_likelihood
#' which is less than 1 will not reach all the way to N * pi_big.
#' Tick marks are drawn on the bars starting from
#' the right end and proceeding at a constant distance that is equal to a single fish
#' with PofZ (top bar) or Scaled Likelihood (bottom bar) of 1.0.  This helps to gauge
#' how far off from 1.0 the individual cells are.  Note that if the sum of all the PofZ's
#' (or of the scaled likelihoods) is less than 1.0, let's say it is f, then the first
#' tick mark will go at the far right side still, but all the other values will be
#' shrunk down so that the sum is f * the max right side value.  This lets us see
#' visually those groups with a cumulative posterior prob across individuals of less than 1.0.
#' @param arrange_scaled_likelihood_ascending If this is TRUE, then the barplot for
#' scaled likelihood has individuals sorted according to their scaled likelihoods.
#' The default is FALSE, in which case the scaled likelihoods are sorted according to
#' the PofZ.  This is better as, I think, as it shows how the interaction of the
#' scaled likelihood and the prior means PofZ is not simply an increasing function
#' of scaled likelihood.
#' If you want to use this, then expand_bars will be set to TRUE.
#' @param dimens named list of relative heights of bars and boxplots. b is relative
#' height of bar, fbox is factor by which the boxplot differs from b, fgap is the
#' fraction of b that is the gap between boxplot and bar and fnext is the number
#' of b's between the bottom of the bottom bar and the top of the next boxplot.
#' @param left_edge_fract the fraction of the sample size to put the lower xlim.
#' Make this larger to keep the N x pi number text from getting cut off.
#' @param Npi_text_size size of text giving the N x pi number
#' @param Npi_left_nudge how much to nudge the N x pi text to the left, away from the
#' start of the bar. Given as a fraction of the sample size.
#' @export
plot_PofZ_vs_pi <- function(
  X,
  mix_coll,
  by_repunit = TRUE,
  burn_in = 100,
  plot_flavor = c("normal", "bars_fully_expanded", "bars_expanded_and_ticked")[1],
  arrange_scaled_likelihood_ascending = FALSE,
  dimens = list(b = 1, fbox = 0.9, fgap = 0.15, fnext = 2),  #
  left_edge_fract = 0.1,
  Npi_text_size = 2,
  Npi_left_nudge = 0.005
) {

  if(!(mix_coll %in% X$mixing_proportions$mixture_collection)) {
    stop("mix_coll ", mix_coll, " not found in X!")
  }
  stopifnot(plot_flavor %in% c("normal", "bars_fully_expanded", "bars_expanded_and_ticked"))


  # capture the input parameters so we can return them as well
  params_list <- as.list(environment())
  # remove the function name and the data set
  params_list$plot_PofZ_vs_pi <- NULL
  params_list$X <- NULL

  # Individual quantity manipulations:
  # 1. Get tibble and filter to mixture collection
  # 2. Calculate scaled likelihoods
  PofZ <- X$indiv_posteriors %>%
    filter(mixture_collection == mix_coll) %>%
    group_by(indiv) %>%
    mutate(
      scaled_likelihood = exp( log_likelihood - max(log_likelihood)),
      scaled_likelihood = scaled_likelihood / sum(scaled_likelihood),
      .after = PofZ
    ) %>%
    ungroup()

  # group by repunit, if indicated
  if(by_repunit == TRUE) {
    PofZ <- PofZ %>%
      group_by(indiv, repunit) %>%
      summarise(
        PofZ = sum(PofZ),
        scaled_likelihood = sum(scaled_likelihood)
      ) %>%
      ungroup() %>%
      rename(group = repunit)
  } else {
    PofZ <- PofZ %>%
      rename(group = collection)
  }



  # Get the mixing proportion posterior means
  Pi_pm <- X$mixing_proportions %>%
    filter(mixture_collection == mix_coll)

  if(by_repunit == TRUE) {
    Pi_pm <- Pi_pm %>%
      group_by(repunit) %>%
      summarise(
        pi = sum(pi)
      ) %>%
      ungroup() %>%
      rename(group = repunit)
  }  else {
    Pi_pm <- Pi_pm %>%
      rename(group = collection)
  }


  # get the mixing proportion traces
  Pi_trace <- X$mix_prop_traces %>%
    filter(mixture_collection == mix_coll, sweep >= burn_in)

  if(by_repunit == TRUE) {
    Pi_trace <- Pi_trace %>%
      group_by(sweep, repunit) %>%
      summarise(
        pi = sum(pi)
      ) %>%
      ungroup() %>%
      rename(group = repunit)
  } else {
    Pi_trace <- Pi_trace %>%
      rename(group = collection)
  }

  # Deal with the vertical dimensions of the bars and things
  D <- dimens$b * (dimens$fbox + dimens$fgap + dimens$fnext + 2)
  bstar <- dimens$b / D
  fbox <- dimens$fbox
  fgap <- dimens$fgap
  fnext <- dimens$fnext

  box_nudge <- bstar * (1 + fgap + 0.5 * fbox)
  top_bar_nudge <- bstar / 2
  bottom_bar_nudge <- -bstar / 2


  # get the order of the groups:
  group_ord <- Pi_pm %>% arrange(pi) %>% pull(group)

  # sample size:
  N <- n_distinct(PofZ$indiv)

  # Prepare orders of groups, etc
  PofZ_ready <- PofZ %>%
    mutate(group_f = factor(group, levels = group_ord)) %>%
    arrange(group_f, PofZ) %>%
    group_by(group_f) %>%
    mutate(spot = n():1, PofZfill = PofZ) %>%
    ungroup() %>%
    mutate(scaled_likelihood_fill = scaled_likelihood)

  ScLik_ready <- PofZ_ready

  Pi_pm_ready <- Pi_pm %>%
    mutate(group_f = factor(group, levels = group_ord)) %>%
    mutate(pipm_times_n = sprintf("%.3f   ", pi * N))  # I am adding a space here to make it easy to fit left of the bar

  Pi_trace_ready <- Pi_trace %>%
    mutate(group_f = factor(group, levels = group_ord))


  # re-arrange scaled likelihoods if called for
  if(arrange_scaled_likelihood_ascending == TRUE) {
    ScLik_ready <- ScLik_ready %>%
      arrange(group_f, scaled_likelihood) %>%
      group_by(group_f) %>%
      mutate(spot = n():1)

  }

  # expand the bars, if indicated
  if(plot_flavor == "bars_fully_expanded") {

    # then we have to figure out the scaling factors for the scaled likelihoods
    # and the posterior probabilities. We want both bars to reach to the same
    # extent across all groups
    scale_factors <- PofZ_ready %>%
      group_by(group) %>%
      summarise(
        sumNPofZ = sum(PofZ),
        sumNScLik = sum(scaled_likelihood)
      ) %>%
      mutate(
        Scale_Factor_PofZ = max(sumNPofZ) / sumNPofZ,
        Scale_Factor_ScLik = max(sumNPofZ) / sumNScLik
      ) %>%
      select(group, Scale_Factor_PofZ, Scale_Factor_ScLik)

    # then we can apply that Scale Factor to each of them
    PofZ_ready <- PofZ_ready %>%
      left_join(scale_factors, by = join_by(group)) %>%
      mutate(
        PofZ = Scale_Factor_PofZ * PofZ,
        scaled_likelihood = Scale_Factor_ScLik * scaled_likelihood
      )
    ScLik_ready <- ScLik_ready %>%
      left_join(scale_factors, by = join_by(group)) %>%
      mutate(
        scaled_likelihood = Scale_Factor_ScLik * scaled_likelihood
      )
  } else if(plot_flavor == "bars_expanded_and_ticked") {
    # now we do something slightly differently for the bars_expanded_and_ticked case.
    scale_factors <- PofZ_ready %>%
      group_by(group) %>%
      summarise(
        sumNPofZ = sum(PofZ),
        sumNScLik = sum(scaled_likelihood)
      ) %>%
      mutate(
        Scale_Factor_PofZ = ifelse(sumNPofZ > 1, max(sumNPofZ) / sumNPofZ, max(sumNPofZ)),
        Scale_Factor_ScLik = ifelse(sumNScLik > 1, max(sumNPofZ) / sumNScLik, max(sumNPofZ))
      ) %>%
      select(group, Scale_Factor_PofZ, Scale_Factor_ScLik)
      # now, get the left-right positions of the ticks
  }


  g <- ggplot() +
    geom_boxplot(
      data = Pi_trace_ready,
      aes(x = group_f, y = pi * N),
      width = (dimens$fbox * dimens$b) / D,
      linewidth = 0.2,
      outlier.stroke = 0.25,
      outlier.size = 0.1,
      position = position_nudge(x = box_nudge)
    ) +
    geom_col(
      data = PofZ_ready,
      mapping = aes(x = group_f, y = PofZ, fill = PofZfill, group = spot),
      width = dimens$b / D,
      colour = "black",
      linewidth = 0.02,
      position = ggpp::position_stacknudge(x = top_bar_nudge)
    ) +
    geom_col(
      data = ScLik_ready,
      mapping = aes(x = group_f, y = scaled_likelihood, fill = scaled_likelihood_fill, group = spot),
      width = dimens$b / D,
      colour = "black",
      linewidth = 0.02,
      position = ggpp::position_stacknudge(x = bottom_bar_nudge)
    ) +
    geom_text(
      data = Pi_pm_ready,
      mapping = aes(x = group_f, label = pipm_times_n),
      y = 0,
      hjust = 1,
      vjust = 0.5,
      size = Npi_text_size
    ) +
    scale_fill_viridis_c() +
    ylim(-left_edge_fract * N, NA) +
    coord_flip() +
    theme_bw()



  # We return the plot, but we also want to return the data that went into
  # making the plot, (in case someone wants to mess with it further)
  # and also the parameters
  list(
    plot = g,
    params = params_list,
    underlying_data = list(
      Pi_pm_ready = Pi_pm_ready,
      Pi_trace_ready = Pi_trace_ready,
      PofZ_ready = PofZ_ready,
      ScLik_ready = ScLik_ready
    )
  )

}

