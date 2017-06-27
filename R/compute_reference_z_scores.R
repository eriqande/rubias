

#' compute empirical Z-scores for reference population individuals
#'
#' For each collection, we compute the mean and sd of the log-likelihoods
#' of assigning correctly to the collection.  We then compute the Z-values
#' for each individual.
#' @param SA The tibble that is output by \code{\link{self_assign}}.
#' @param returnCollectionSummary if this is TRUE then this just returns the variance and the L_bar
#' for each collection so that it can be used to compute z-scores for mixture fish
#' @return A tibble with a few more columns
#' @export
compute_reference_z_scores <- function(SA, returnCollectionSummary = FALSE) {
  tmp <- SA %>%
    dplyr::filter(collection == inferred_collection) %>%
    dplyr::group_by(collection) %>%
    dplyr::mutate(mean = mean(log_likelihood),
                  sd = sd(log_likelihood),
                  naive_z = (log_likelihood - mean) / sd,
                  L_bar = sum(log_likelihood) / sum(n_non_miss_loci),
                  SS = sum(ssq_logl),
                  nn = sum(n_non_miss_loci)) %>%
    dplyr::mutate(var1 = (1 / (nn - 1)) * (SS - nn * (L_bar ^ 2)),
                  z_score = (log_likelihood - (n_non_miss_loci * L_bar)) / sqrt(n_non_miss_loci * var1))

  if (returnCollectionSummary == TRUE) {
    ret <- tmp %>%
      dplyr::summarise(L_bar = mean(L_bar),
                var1 = mean(var1))  %>%
      dplyr::ungroup() %>%
      dplyr::select(collection, L_bar, var1)
  } else {
    ret <- tmp %>%
      dplyr::select(-mean, -sd, -L_bar, -SS, -nn, -var1)
  }
  ret
}
