

#' A helper function to check that the input data frame is OK
#'
#' Just checks to make sure that column types are correct
#' @param D the data frame
#' @param gen_start_col  the column in which the genetic data starts
#' @param type For writing errors, supply "mixture" or "reference" as appropriate.
#' @keywords internal
check_refmix <- function(D, gen_start_col, type = "reference") {

  # first check to make sure that the repunit, collection, and indiv columns are present
  if(!("repunit") %in% names(D)) stop("Missing column \"repunit\" in", type)
  if(!("collection") %in% names(D)) stop("Missing column \"collection\" in", type)
  if(!("indiv") %in% names(D)) stop("Missing column \"indiv\" in", type)

  # now check to see if any of those are not character vectors
  if(!is.character(D$repunit)) stop("Column \"repunit\" must be a character vector.  It is not in ", type, " data frame")
  if(!is.character(D$collection)) stop("Column \"collection\" must be a character vector.  It is not in ", type, " data frame")
  if(!is.character(D$indiv)) stop("Column \"indiv\" must be a character vector.  It is not in ", type, " data frame")

  # now, check to make sure that all the locus columns are character or integer:
  tmp <- D[, -(1:(gen_start_col - 1))]
  char_or_int <- sapply(tmp, is.character) | sapply(tmp, is.integer)
  if(any(!char_or_int)) {
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
#' @param coll_levs a vector of levels for the collection (or NULL if none)
#' @param repu_levs a vector of levels for the repunit (or NULL if none)
#' @keywords internal
tidy_mcmc_coll_rep_stuff <- function(field, p, pname, car_tib, coll_levs = NULL, repu_levs = NULL) {
  ret <- tibble::tibble(collection = car_tib$collection, value = field[[p]]) %>%
    dplyr::left_join(car_tib, ., by = "collection")

  # change the name
  names(ret)[names(ret) == "value"] <- pname

  if (!is.null(coll_levs)) {
    ret$collection = factor(ret$collection, levels = coll_levs)
  }
  if (!is.null(repu_levs)) {
    ret$repunit = factor(ret$repunit, levels = repu_levs)

  }
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
#' @param coll_levs a vector of levels for the collection (or NULL if none)
#' @param repu_levs a vector of levels for the repunit (or NULL if none)
#' @keywords internal
tidy_mcmc_pofz <- function(input, pname, car_tib, mix_indiv_tib, coll_levs = NULL, repu_levs = NULL) {
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

  # deal with the factor levels if present
  if (!is.null(coll_levs)) {
    ret$collection = factor(ret$collection, levels = coll_levs)
  }
  if (!is.null(repu_levs)) {
    ret$repunit = factor(ret$repunit, levels = repu_levs)

  }
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
#' @param coll_levs a vector of levels for the collection (or NULL if none)
#' @param repu_levs a vector of levels for the repunit (or NULL if none)
#' @param interval the thinning interval that was used
#' @keywords internal
tidy_pi_traces <- function(input, pname, car_tib, coll_levs, repu_levs, interval) {
  ret <- lapply(input, function(x) tibble::tibble(collection = car_tib$collection,
                                            pi = x)) %>%
    dplyr::bind_rows(.id = "sweep") %>%
    dplyr::mutate(sweep = as.integer(sweep) * as.integer(interval)) %>%
    dplyr::left_join(., car_tib, by = "collection") %>%
    dplyr::select(sweep, repunit, collection, pi)

  names(ret)[names(ret) == "pi"] <- pname

  if (!is.null(coll_levs)) {
    ret$collection = factor(ret$collection, levels = coll_levs)
  }
  if (!is.null(repu_levs)) {
    ret$repunit = factor(ret$repunit, levels = repu_levs)
  }
  ret
}




#' Generate a samples for a mixture.
#'
#' Creates random reporting unit (rho) and collection (omega) proportions, and a
#' \code{sim_colls} vector for simulation of individual genotypes, based on the methods
#' used in Hasselman \emph{et al.} (2015)
#'
#' @param RU_starts a vector delineating the reporting units in \code{RU_vec};
#' generated by \code{tcf2param_list}
#' @param RU_vec a vector of collection indices, grouped by reporting unit;
#' generated by \code{tcf2param_list}
#' @param size The number of individuals desired in the mixture sample.  Default = 100.
#' @param alpha_repunit The dirichlet parameter for simulating the proportions of reporting units. Default = 1.5
#' @param alpha_collection The dirichlet parameter for simulating proportions of collections within reporting units. Default = 1.5
#'
#' @return a list with three elements.
#' The first two are a rho vector and an omega vector, respectively. The third is a vector of origins for
#' simulated individuals, sampled from the collections with probabilities = omega
#' @export
#' @keywords internal
simulate_random_samples <- function(RU_starts, RU_vec, size = 100, alpha_repunit = 1.5, alpha_collection = 1.5) {
  rho <- gtools::rdirichlet(1, rep(alpha_repunit, length(RU_starts) - 1))
  omega <- numeric(length(RU_vec))
  omega <- sapply(1:(length(RU_starts)-1), function(x) {
    omega[RU_vec[(RU_starts[x] + 1):RU_starts[x+1]]] <- gtools::rdirichlet(1, rep(alpha_collection, RU_starts[x + 1] - RU_starts[x])) * rho[x]
  }) %>%
    cbind() %>%
    unlist()
  sim_coll = sample(RU_vec, size = size, replace = T, prob = omega)

  # finally return it all in a list
  list(rho = rho, omega = omega, sim_coll = sim_coll)
}




