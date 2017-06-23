

#' A helper function to check that the input data frame is OK
#'
#' Just checks to make sure that column types are correct
#' @param D the data frame
#' @param gen_start_col  the column in which the genetic data starts
#' @param type For writing errors, supply "mixture" or "reference" as appropriate.
#' @keywords internal
#' @export
check_refmix <- function(D, gen_start_col, type = "reference") {

  # first check to make sure that the repunit, collection, and indiv columns are present
  if (!("repunit") %in% names(D)) stop("Missing column \"repunit\" in", type)
  if (!("collection") %in% names(D)) stop("Missing column \"collection\" in", type)
  if (!("indiv") %in% names(D)) stop("Missing column \"indiv\" in", type)

  # now check to see if any of those are not character vectors
  if (!is.character(D$repunit)) stop("Column \"repunit\" must be a character vector.  It is not in ", type, " data frame")
  if (!is.character(D$collection)) stop("Column \"collection\" must be a character vector.  It is not in ", type, " data frame")
  if (!is.character(D$indiv)) stop("Column \"indiv\" must be a character vector.  It is not in ", type, " data frame")

  # now, check to make sure that all the locus columns are character or integer:
  tmp <- D[, -(1:(gen_start_col - 1))]
  char_or_int <- sapply(tmp, is.character) | sapply(tmp, is.integer)
  if (any(!char_or_int)) {
    stop("All locus columns must be of characters or integers.  These in ", type, " are not: ",
         paste(names(char_or_int[!char_or_int]), collapse = ", "))
  }

}


#' A helper function to tidy up the output from the gsi_mcmc functions
#'
#' This makes a tidy data frame of stuff, and also changes things back to
#' factors, if the levels are provided.
#' @param field  The output to tidy (i.e.. out$mean)
#' @param p the name of the parameter whose values you want to extract (like "pi")
#' @param pname the name that you want the parameter to be called in the output
#' @param car_tib  a tibble with repunit and collection in the order they appear in the output
#' @keywords internal
tidy_mcmc_coll_rep_stuff <- function(field, p, pname, car_tib) {
  ret <- tibble::tibble(collection = car_tib$collection, value = field[[p]]) %>%
    dplyr::left_join(car_tib, ., by = "collection")

  # change the name
  names(ret)[names(ret) == "value"] <- pname

  ret
}






#' A helper function to tidy up the PofZ-like output from the gsi_mcmc functions
#'
#' This makes a tidy data frame of stuff, and also changes things back to
#' factors, if the levels are provided.
#' @param input  The output to tidy (i.e.. out$mean$PofZ)
#' @param pname the name that you want the parameter to be called in the output
#' @param car_tib  a tibble with repunit and collection in the order they appear in the output
#' @param mix_indiv_tib  a tibble with the individuals in the order they appear in the output
#' @keywords internal
tidy_mcmc_pofz <- function(input, pname, car_tib, mix_indiv_tib) {
  pofz_mat <- t(input)
  colnames(pofz_mat) <- car_tib$collection

  ret <- dplyr::bind_cols(mix_indiv_tib,
                                tibble::as_tibble(pofz_mat)) %>%
    tidyr::gather(data = ., key = "collection", value = "pofz", -indiv) %>%
    dplyr::left_join(car_tib, by = "collection") %>%
    dplyr::select(indiv, repunit, collection, pofz) %>%
    dplyr::mutate(repunit = factor(repunit, levels = unique(car_tib$repunit)),   # this is heinous uglyness to get populations sorted in the order they came in the data set if they aren't factors to begin with
                  collection = factor(collection, levels = car_tib$collection)) %>%
    dplyr::arrange(indiv, repunit, collection) %>%
    dplyr::mutate(repunit = as.character(repunit),
                  collection = as.character(collection)) %>%
    dplyr::left_join(mix_indiv_tib, ., by = "indiv")

  names(ret)[names(ret) == "pofz"] <- pname


  ret
}


#' a helper function to tidy up the pi-traces that come out of the mcmc functions
#'
#' This makes a tidy data frame of stuff, and also changes things back to
#' factors, if the levels are provided.
#' @param input  The output to tidy (i.e.. out$trace$pi)
#' @param pname the name that you want the parameter to be called in the output
#' @param car_tib  a tibble with repunit and collection in the order they appear in the output
#' @param mix_indiv_tib  a tibble with the individuals in the order they appear in the output
#' @param interval the thinning interval that was used
#' @keywords internal
tidy_pi_traces <- function(input, pname, car_tib, interval) {
  ret <- lapply(input, function(x) tibble::tibble(collection = car_tib$collection,
                                            pi = x)) %>%
    dplyr::bind_rows(.id = "sweep") %>%
    dplyr::mutate(sweep = as.integer(sweep) * as.integer(interval)) %>%
    dplyr::left_join(., car_tib, by = "collection") %>%
    dplyr::select(sweep, repunit, collection, pi)

  names(ret)[names(ret) == "pi"] <- pname

  ret
}





