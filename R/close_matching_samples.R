
#' check for matching (or close to matching) genotypes in a data frame
#'
#' Super simple function that looks at all pairs of fish from the data frame and
#' returns a tibble that includes those which shared a fraction >= than min_frac_non_miss
#' of the genotypes not missing in either fish, and which were matching at a fraction >= min_frac_matching
#' of those non-missing pairs of genotypes.
#' @param D a two-column format genetic dataset, with "repunit", "collection", and "indiv"
#' columns, as well as a "sample_type" column that has entried either of "reference" or "mixture"
#' or both.
#' @param min_frac_non_miss the fraction of loci that the pair must share non missing in order to be reported
#' @param min_frac_matching the fraction of shared non-missing loci that must be shared between the indivdiuals
#' to be reported as a matching pair.
#' @inheritParams assess_reference_loo
#' @return a tibble ...
#' @export
#' @examples
#' # one pair found in the interal alewife data set:
#' close_matching_samples(alewife, 17)
close_matching_samples <- function(D, gen_start_col, min_frac_non_miss = .7, min_frac_matching = .9) {

    # make sure the reference file is OK
    ploidies <- check_refmix(D, gen_start_col, "reference")


    # get the necessary parameters from the reference data
    params <- tcf2param_list(D,
                             gen_start_col,
                             summ = T,
                             alle_freq_prior = list("const_scaled" = 1),
                             ploidies = ploidies)

    # now, we have a function written in Rcpp that does the pairwise comparisons
    # using the information stored in params.  It returns a list of vectors
    # with indices of the matching individuals.
    matches <- tibble::as_tibble(rcpp_close_matchers(params, min_frac_non_miss, min_frac_matching))

    # now we join the individual names
    D1 <- D %>%
      dplyr::ungroup() %>%
      dplyr::select(sample_type, repunit, collection, indiv) %>%
      setNames(paste0(names(.), "_1")) %>%
      dplyr::mutate(idx = 1:n())
    D2 <- D %>%
      dplyr::ungroup() %>%
      dplyr::select(sample_type, repunit, collection, indiv) %>%
      setNames(paste0(names(.), "_2")) %>%
      dplyr::mutate(idx = 1:n())

    # then reorder the columns and return em
    matches %>%
      dplyr::left_join(D1, by = c("indx1" = "idx")) %>%
      dplyr::left_join(D2, by = c("indx2" = "idx")) %>%
      dplyr::select(-indx1, -indx2) %>%
      dplyr::select(num_non_miss, num_match, indiv_1, indiv_2, collection_1, collection_2, dplyr::everything())





}
