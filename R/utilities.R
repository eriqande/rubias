#' Separate a chosen proportion of a reference dataset into a mixture with known population proportions
#'
#' Takes a reference dataset and a set of population proportions, either at the collection
#' or reporting unit level. Randomly samples individuals to satisfy these desired proportions,
#' and splits them into a new "mixture" dataframe.
#'
#'
#' @param D a two-column genetic dataframe with "indiv", "repunit", and "collection" columns
#' @param rhos a vector of the desired reporting unit proportions in the mixture set;
#' if not named, will be assumed to be ordered by order of appearance in the dataset
#' @param omegas the desired collection proportions in the mixture set
#' @param N the total size of the mixture set
#' @param min_remaining the fraction of any collection in the reference dataset which must remain
#' at the end of the draw
#'
#' @return \code{mixture_draw} returns a list of two data frames,
#' "mixture" being the random sample taken, and "reference" being the remaining samples
#'
#' @examples
#' rhos <- as.vector(gtools::rdirichlet(1, table(alewife$repunit)))
#' cross_val <- mixture_draw(D = alewife, rhos = rhos, N = 100, min_remaining = .005)
#'
#' @export
mixture_draw <- function(D, rhos = NULL, omegas = NULL, N, min_remaining = 1/(length(omegas) * 10)) {
  if(!is.null(rhos) && !is.null(omegas)) stop("Cannot specify proportions of both rho and omega")
  repidxs <- D %>%
    dplyr::select(repunit, collection) %>%
    dplyr::group_by(repunit, collection) %>%
    dplyr::tally()
  coll_props <- repidxs$n /sum(repidxs$n)

  if(sum(repidxs$n) < N) stop("Cannot take mixture sample of size greater than reference dataset")
  if(any(coll_props < min_remaining)) stop("one or more collections starting with proportion below the specified minimum")


  # deterministic sampling (rho is exactly the desired value)
  if(!is.null(rhos)) {
    if(!identical(all.equal(sum(rhos),1), TRUE)) stop("Desired proportions must sum to 1")

    if(is.null(names(rhos))) names(rhos) <- levels(as.factor(D$repunit))
    ru_ns <- table(D$repunit)
    samp_sizes <- as.vector(round2(rhos * N, 0))
    names(samp_sizes) <- names(rhos)
    draw <- lapply(names(rhos), function(repunit) {
      samp <- dplyr::sample_n(D[D$repunit == repunit, ], samp_sizes[repunit], replace = FALSE)
    }) %>%
      do.call("rbind", .)
    new_D <- dplyr::setdiff(D, draw)
    new_repidxs <- new_D %>%
      dplyr::select(repunit, collection) %>%
      dplyr::group_by(repunit, collection) %>%
      dplyr::tally()
    coll_props <- new_repidxs$n /sum(new_repidxs$n)
  }

  else if(!is.null(omegas)) {
    if(!identical(all.equal(sum(omegas),1), TRUE)) stop("Desired proportions must sum to 1")

    if(is.null(names(omegas))) names(omegas) <- levels(D$collection)
    samp_sizes <- round2(omegas * N, 0)
    draw <- lapply(levels(repidxs$collection), function(collection) {
      samp <- dplyr::sample_n(D[D$collection == collection, ], samp_sizes[collection], replace = FALSE)
    }) %>%
      do.call("rbind", .)
    new_D <- dplyr::setdiff(D, draw)
    new_repidxs <- new_D %>%
      dplyr::select(repunit, collection) %>%
      dplyr::group_by(repunit, collection) %>%
      dplyr::tally()
    coll_props <- new_repidxs$n /sum(new_repidxs$n)
  }

  draw$sample_type <- "mixture"
  new_D$sample_type <- "reference"
  out <- list("mixture" = draw, "reference" = new_D)

}


#' Round a given number, with 5 always rounded up
#'
#' Given a number and a digit to round to, returns the rounded number,
#' with 5 always rounded upwards.
#'
#' This function differs from \code{round}, which rounds 5 "towards the
#' even number". Rounding 5s up leads to bias when positive and negative numbers
#' are expected, but can be desired in some cases.
#'
#' @param x the data to be rounded
#' @param n the number of digits to round to
#'
#' @export
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}
