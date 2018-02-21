

#' for each individual count the number or loci missing and non_missing
#'
#' Takes a rubias genetic data frame that must have column "indiv". Haploids have the
#' second column at each locus totally missing.  Diploids with missing data will have both gene copies
#' missing.
#' @param D the data frame
#' @param gen_start_col the column in which the genetic data starts
#' @return returns a data frame with indiv (as characters), n_non_miss_loci, n_miss_loci (as numeric) and missing_loci
#' (as a list-column of named integer vectors)
#' @keywords internal
#' @export
count_missing_data <- function(D, gen_start_col) {

  D2 <- D[, seq(gen_start_col, ncol(D), by = 2)]  # just take the first column for each locus

  miss <- as.integer(rowSums(is.na(D2)))
  nonmiss <- as.integer(rowSums(!is.na(D2)))

  # finally, let's also store the pattern of missing data in a list column
  miss_pattern_list <- apply(D2, 1, function(x) which(is.na(x)))
  # was having problems when no missing data is present; output to missing pattern list
  # is a single empty integer, which is not accepted for tibble creation. Introduced
  # the following fix:
  if(length(miss_pattern_list) == 0) miss_pattern_list <- list(integer())[rep(1,nrow(D))]


  tibble::tibble(indiv = as.character(D$indiv), n_non_miss_loci = nonmiss, n_miss_loci = miss, missing_loci = miss_pattern_list)

}
