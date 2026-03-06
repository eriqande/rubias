#' Calculates the posterior probability that the sample contains at least one fish from each collection and repunit
#'
#' This takes the output of infer_mixture() and returns two tibbles, one for collections
#' and the other for repunits. Note that for this to work properly, infer_mixture() must
#' have been run with a total_catch_tib, and also it is imperative that the
#' variable_prob_is_catch parameter to infer_mixture() was at its default value of
#' FALSE.
#' @param X The return object from infer_mixture
#' @param mix_coll  character name of the mixture_collection within X
#' you want to summarize the results for.
#' @param burn_in  discards the first burn_in sweeps when handling the pi_traces
#' (for the boxplots of the pi)
#' @export
prob_at_least_one_in_sample <- function(
    X,
    mix_coll,
    burn_in = 100
) {

  if(!(mix_coll %in% X$mixing_proportions$mixture_collection)) {
    stop("mix_coll ", mix_coll, " not found in X!")
  }

  if(!("allocation_count_traces" %in% names(X) )) {
    stop("X must come from infer_mixture() run with a total_catch_tib to use this function.")
  }

  # here we get the by-collection version
  by_coll <- X$allocation_count_traces %>%
    filter(mixture_collection == mix_coll, sweep >= burn_in) %>%
    group_by(collection) %>%
    summarise(fract_non_zero = mean(CA > 0)) %>%
    ungroup() %>%
    arrange(desc(fract_non_zero))

  # here we get the by-repunit version
  by_repu <- X$allocation_count_traces %>%
    filter(mixture_collection == mix_coll, sweep >= burn_in) %>%
    group_by(sweep, repunit) %>%
    summarise(CA_repu = sum(CA)) %>%
    ungroup() %>%
    group_by(repunit) %>%
    summarise(fract_non_zero = mean(CA_repu > 0)) %>%
    ungroup() %>%
    arrange(desc(fract_non_zero))

  list(
    by_repunit = by_repu,
    by_collection = by_coll
  )
}
