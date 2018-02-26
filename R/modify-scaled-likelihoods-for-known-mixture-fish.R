

#' for individuals of known origin in the mixture, put all their weight on their known collection
#'
#' This is used internally.
#' @param SL the matrix of scaled likelihoods.
#' @param KC a character vector of collections that the individuals belong to (or NA if you
#' don't know where they come from).  If this is NULL, then SL just gets returned untouched.
#' @param CFL the levels of the collections factor (which is within clean$short)
#' @export
#' @keywords internal
modify_scaled_likelihoods_for_known_mixture_fish <- function(SL, KC, CFL) {
  # if there is no known collection column, just return SL untouched
  if(is.null(KC)) {
    return(SL)
  }
  ret <- SL
  # otherwise, do what needs to be done
  ipops <- as.integer(factor(KC, levels = CFL))
  for(i in 1:ncol(SL)) {
    if(!is.na(ipops[i])) {
      ret[, i] <- 1e-100  # putting this in here in case we end up taking logs of it or something. (1 over a googol is close to 0...)
      ret[ipops[i], i] <- 1.0
    }
  }
  ret
}
