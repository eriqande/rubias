#' Make a plot comparing PofZs and Pi values for all collections (or repunits)
#'
#' @param X The return object from infer_mixture
#' @param mix_coll  character name of the mixture_collection within X
#' you want to plot results for.
#' @param by_repunit logical.  If true, do by reporting units rather than collections.
#' @param burn_in  discards the first burn_in sweeps when handling the pi_traces
#' (for the boxplots of the pi)
#' @export
plot_PofZ_vs_pi <- function(
  X,
  mix_coll,
  by_repunit = TRUE,
  burn_in = 100,
  dimens = list(b = 1, fbox = 0.9, fgap = 0.15, fnext = 2)  # b is relative height of bar, fbox is factor by
     # which the boxplot differs from b, fgap is the fraction of b that is the gap between boxplot and bar
     # and fnext is the number of b's between the bottom of the bottom bar and the top of the next boxplot.
) {

  if(!(mix_coll %in% X$mixing_proportions$mixture_collection)) {
    stop("mix_coll ", mix_coll, " not found in X!")
  }

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

  full_nudge <- 1 / 3

  first_bar_center <- full_nudge - (0.5 * (dimens$b  * dimens$fbox) + dimens$b * dimens$fgap + 0.5 * dimens$b) / D
  first_bar_top <- first_bar_center + (0.5 * dimens$b) / D
  second_bar_center <- first_bar_center - dimens$b / D

  # get the order of the groups:
  group_ord <- Pi_pm %>% arrange(pi) %>% pull(group)

  # sample size:
  N <- n_distinct(PofZ$indiv)

  # testing geom_col
  PofZ_ready <- PofZ %>%
    mutate(group_f = factor(group, levels = group_ord)) %>%
    arrange(group_f, PofZ) %>%
    group_by(group_f) %>%
    mutate(spot = n():1, PofZfill = PofZ) %>%
    ungroup() %>%
    mutate(scaled_likelihood_fill = scaled_likelihood)

  Pi_trace_ready <- Pi_trace %>%
    mutate(group_f = factor(group, levels = group_ord))

  ggplot() +
    geom_boxplot(
      data = Pi_trace_ready,
      aes(x = group_f, y = pi * N),
      width = (dimens$fbox * dimens$b) / D,
      linewidth = 0.2,
      outlier.stroke = 0.25,
      outlier.size = 0.1,
      position = position_nudge(x = full_nudge)
    ) +
    geom_col(
      data = PofZ_ready,
      mapping = aes(x = group_f, y = PofZ, fill = PofZfill, group = spot),
      width = dimens$b / D,
      colour = "black",
      linewidth = 0.02,
      position = ggpp::position_stacknudge(x = first_bar_center)
    ) +
    geom_col(
      data = PofZ_ready,
      mapping = aes(x = group_f, y = scaled_likelihood, fill = scaled_likelihood_fill, group = spot),
      width = dimens$b / D,
      colour = "black",
      linewidth = 0.02,
      position = ggpp::position_stacknudge(x = second_bar_center)
    ) +
    scale_fill_viridis_c() +
    coord_flip() +
    theme_bw()


    geom_col(position = "stack", colour = "black", linewidth = 0.02) +

    coord_flip()

}

