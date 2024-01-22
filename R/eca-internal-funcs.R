
#' infer the ploidy from the pattern of NAs in the columns of data
#'
#' This is strictly internal
#' @param tmp a data frame with 2 * L columns (two for each locus)
#' @keywords internal
get_ploidy_from_frame <- function(tmp, type) {
  ploidy <- rep(2, ncol(tmp) / 2)  # initialize to diploid
  locus_names <- names(tmp)[c(TRUE, FALSE)]
  gc_mism_error <- FALSE
  j <- 0
  for (i in seq(1, ncol(tmp), by = 2)) {
    a <- tmp[,i]
    b <- tmp[,i+1]
    j <- j + 1

    looksHaploid <- all(is.na(b))
    gc_mism <- any(xor(is.na(a), is.na(b)))  # returns true if one gene copy is missing and not the other for any locus

    if (looksHaploid == TRUE) {
      if(any(!is.na(a))) {
        ploidy[j] <- 1
        message("Scoring locus ", names(tmp)[i], " as haploid")
      } else {
        if(all(is.na(a))) {
          ploidy[j] <- 0
          message("All gene copies missing at locus ", names(tmp)[i], " in ", type, " data. Ploidy indeterminate in that data frame.")
        }
      }
    } else {
      if(gc_mism == TRUE) {
        message("Error in input.  At diploid loci, either both or neither gene copies must be missing. Offending locus = ", names(tmp)[i], "\n\tNote! This might indicate that the gen_start_col is incorrect.")
        gc_mism_error <- TRUE
      }
    }
  }
  if(gc_mism_error == TRUE) stop("Bailing out due to single gene copies being missing data at non-haploid loci.")

  names(ploidy) <- locus_names
  ploidy
}


#' A helper function to check that the input data frame is OK
#'
#' Just checks to make sure that column types are correct.
#'
#' It also checks the patterns of missing data, and from that infers whether
#' markers are haploid or diploid.
#' @param D the data frame
#' @param gen_start_col  the column in which the genetic data starts
#' @param type For writing errors, supply "mixture" or "reference" as appropriate.
#' @keywords internal
#' @export
check_refmix <- function(D, gen_start_col, type = "reference") {

  # first make sure that the data frame is not grouped
  if (dplyr::is_grouped_df(D) == TRUE) {
    stop("D is a grouped data frame.  Please ungroup it before using it in rubias.")
  }

  # first check to make sure that the repunit, collection, and indiv columns are present
  if (!("repunit") %in% names(D)) stop("Missing column \"repunit\" in", type)
  if (!("collection") %in% names(D)) stop("Missing column \"collection\" in", type)
  if (!("indiv") %in% names(D)) stop("Missing column \"indiv\" in", type)

  # check to make sure that if any fish has sample_type is "mixture" then it has NA for repunit
  mixture_nonNAs <- D %>%
    dplyr::select(sample_type, repunit) %>%
    dplyr::filter(sample_type == "mixture" & !is.na(repunit))


  # now check to see if any of those are not character vectors.  Mixture samples that have NA for their
  # repunits  might be logicals, so we only freak out if they aren't all NAs and are still not characters
  if (!all(is.na(D$repunit)) && !is.character(D$repunit)) stop("Column \"repunit\" must be a character vector.  It is not in ", type, " data frame")
  if (!is.character(D$collection)) stop("Column \"collection\" must be a character vector.  It is not in ", type, " data frame")
  if (!is.character(D$indiv)) stop("Column \"indiv\" must be a character vector.  It is not in ", type, " data frame")


  if(nrow(mixture_nonNAs) > 0) {
    stop("Error. All fish of sample_type == \"mixture\" must have repunit == NA. ",
         nrow(mixture_nonNAs),
         " fish have sample_type == \"mixture\" but repunit is: ",
         paste(unique(mixture_nonNAs$repunit), collapse = ", "))
  }

  # now, check to make sure that all the locus columns are character or integer (or logical, as the case
  # might be if they are all missing to denote haploids).
  tmp <- D[, -(1:(gen_start_col - 1))]
  char_or_int <- sapply(tmp, is.character) | sapply(tmp, is.integer) | sapply(tmp, is.logical)

  if (any(!char_or_int)) {
    stop("All locus columns must be of characters or integers.  These in ", type, " are not: ",
         paste(names(char_or_int[!char_or_int]), collapse = ", "))
  }

  # now cycle over the loci and check the pattern of missing data.  Any individual
  # with missing data must be missing at both gene copies, unless it is a haploid
  # marker, in which case it must be missing at the second gene copy in everyone.
  ploidy <- get_ploidy_from_frame(tmp, type = type)


  # check also to make sure that indiv IDs are unique
  dupies <- tibble::tibble(indiv = D$indiv) %>%
    dplyr::count(indiv) %>%
    dplyr::filter(n > 1)
  if (nrow(dupies) > 0) {
    err_string <- paste(sprintf("%s (%d)", dupies$indiv, dupies$n), collapse = ", ")
    stop("The following indiv IDs are repeated.  Number of times in parentheses: ", err_string)
  }

  # now, check to make sure that no collection is associated with more than one repunit
  msc <- D[, c("repunit", "collection")] %>%
    dplyr::count(repunit, collection) %>%
    dplyr::filter(n > 0) %>%
    dplyr::group_by(collection) %>%
    dplyr::mutate(times_seen = n()) %>%
    dplyr::filter(times_seen > 1) %>%
    dplyr::arrange(collection)

  if (nrow(msc) > 0) {
    err_str <- paste(paste(msc$n, "indivs in collection", msc$collection, "in repunit", msc$repunit ), collapse = ",\n\t")
    stop("Each collection must belong to no more than one repunit.  Offenders in ", type, ":\n\t", err_str)
  }

  # at the end, we return a vector of ploidies (1 or 2) for the loci
  ploidy
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
#' @export
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
#' @export
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
#' @param interval the thinning interval that was used
#' @keywords internal
#' @export
tidy_pi_traces <- function(input, pname, car_tib, interval) {


  ret <- tibble::tibble(
    sweep = rep(0:(length(input) - 1), each = length(input[[1]])),
    collection = rep(car_tib$collection, length(input)),
    pi = unlist(input)
  ) %>%
    dplyr::mutate(sweep = as.integer(sweep) * as.integer(interval)) %>%
    dplyr::left_join(., car_tib, by = "collection") %>%
    dplyr::select(sweep, repunit, collection, pi)

  names(ret)[names(ret) == "pi"] <- pname

  ret
}





