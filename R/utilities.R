#' Separate a chosen proportion of a reference dataset into a mixture with known population proportions
#'
#' @param D a two-column genetic dataframe with "indiv", "repunit", and "collection" columns
#' @param rhos the desired reporting unit proportions in the mixture set;
#' if not named, will be assumed to be ordered by order of appearance in the dataset
#' @param omegas the desired collection proportions in the mixture set
#' @param N the total size of the mixture set
#' @param min_remaining the fraction of any collection in the reference dataset which must remain
#' at the end of the draw
#'
#' @return a list of two data frames, "mixture" being the random sample taken, and "reference" being the remaining samples
#'
#' D <- alewife
#' N <- 100
#' rhos <- as.vector(gtools::rdirichlet(1, table(alewife$repunit)))
#' min_remaining <- .005
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

    if(is.null(names(rhos))) names(rhos) <- levels(D$repunit)
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
    if(any(coll_props < min_remaining)) stop("minimum remaining violated")
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
    if(any(coll_props < min_remaining)) stop("desired sample cannot be taken without violating min_remaining")
  }

  draw$sample_type <- "mixture"
  new_D$sample_type <- "reference"
  out <- list("mixture" = draw, "reference" = new_D)

}


#' Rounding 5 up, rather than rounding to the nearest even number
#'
#' Sometimes desired, especially when the numbers are all positive,
#' and not possible with R (stolen from the Statistically Significant blog)
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
