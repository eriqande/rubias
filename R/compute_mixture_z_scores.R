


#' compute empirical Z-scores for mixture sample population individuals
#'
#' To run this, you must supply the z-scores and logls of the baseline fish, too.
#' @param M The tibble that is output by \code{\link{infer_mixture}}.
#' @param SA The tibble that is returned by \code{\link{self_assign}} using the same reference.
#' @return A tibble like with a few more columns
#' @export
compute_mixture_z_scores <- function(SA) {
  # this is incomplete.  I think that I need to return some more things from the reference z function first

  # for now just return nonsense:
  1
}
