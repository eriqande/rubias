#' Create a vector of pi Dirichlet priors with specified values for one or more collections
#'
#' This handles a case in which the user provides a data frame for \code{pi_prior}. The
#' data frame lists desired Dirichlet parameter priors for at least one reference collection,
#' and/or a default value for all unspecified collections.
#'
#' Input checking is currently done in the early stages of \code{infer_mixture} in order to
#' throw errors before long processing times, and avoid re-checking during \code{bootstrap_rho}.
#'
#' @param P  A data frame of one or more desired pi prior parameters. One column, "collection",
#' is a character vector, with valid values including the names of any reference collections,
#' or the special value "DEFAULT_PI". The second column, "pi_param" is the prior value to be
#' used for each collection.
#' @param C a character vector of all unique reference collection names
#' @keywords internal
#' @export
custom_pi_prior <- function(P, C) {
  if(!("DEFAULT_PI") %in% P$collection) {
    default <- 1/length(C)
  } else default <- P$pi_param[P$collection == "DEFAULT_PI"]
  ret <- dplyr::left_join(C, P) %>%
    tidyr::replace_na(list(pi_param = default))
  ret$pi_param
}

