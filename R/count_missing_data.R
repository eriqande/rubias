

#' for each individual count the number or loci missing and non_missing
#'
#' Takes a rubias genetic data frame that must have column "indiv".  Currently
#' hard-wired for diploids.
#' @param D the data frame
#' @param gen_start_col the column in which the genetic data starts
#' @return returns a data frame with indiv (as characters), n_non_miss_loci, n_miss_loci (as numeric) and missing_loci
#' (as a list-column of named integer vectors)
#' @keywords internal
#' @export
count_missing_data <- function(D, gen_start_col) {

  DM <- as.matrix(D[, gen_start_col:ncol(D)])
  rownames(DM) <- D$indiv

  miss <- rowSums(is.na(DM)) / 2
  nonmiss <- rowSums(!is.na(DM)) / 2

  # finally, let's also store the pattern of missing data in a list column
  D2 <- D[, seq(1, ncol(D), by = 2)]  # just take the first column for each locus
  miss_pattern_list <- apply(D2, 1, function(x) which(is.na(x)))


  tibble::tibble(indiv = as.character(D$indiv), n_non_miss_loci = nonmiss, n_miss_loci = miss, missing_loci = miss_pattern_list)

}
