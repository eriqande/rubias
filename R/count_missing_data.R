

#' for each individual count the number or loci missing and non_missing
#'
#' Takes a rubias genetic data frame that must have column "indiv".  Currently
#' hard-wired for diploids.
#' @param D the data frame
#' @param gen_start_col the column in which the genetic data starts
#' @return returns a data frame with indiv (as characters), n_non_miss_loci, n_miss_loci (as numeric)
#' @keywords internal
count_missing_data <- function(D, gen_start_col) {

  DM <- as.matrix(D[, gen_start_col:ncol(D)])
  rownames(DM) <- D$indiv

  miss <- rowSums(is.na(DM)) / 2
  nonmiss <- rowSums(!is.na(DM)) / 2

  tibble::tibble(indiv = as.character(D$indiv), n_non_miss_loci = nonmiss, n_miss_loci = miss)

}
