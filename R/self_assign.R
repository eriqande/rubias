#' Do leave-one-out self-assignment of individuals in a reference baseline
#'
#' Returns a tidy data frame
#' @inheritParams assess_reference_loo
#' @return a tibble ...
#' @export
#' @examples
#' ale_sa <- self_assign(alewife, 17)
self_assign <- function(reference, gen_start_col) {

  # make sure the reference file is OK
  check_refmix(reference, gen_start_col, "reference")

  # get the necessary parameters from the reference data
  params <- tcf2param_list(reference, gen_start_col, summ = T)

  # get the log-likelihoods
  logl <- t(geno_logL(par_list = params))
  # and get the sum-of-squares over loci of the logls (for z-score calculation)
  logl_ssq <- t(geno_logL_ssq(par_list = params))

  # put the collection names at the top of them. To do this, we put RU_vec into sorted
  # order and then grab the names off it
  colnames(logl) <- names(sort(params$RU_vec))
  colnames(logl_ssq) <- names(sort(params$RU_vec))

  # then make a tibble of the logls and put the meta data (indiv, collection, repuunit) from
  # "reference" back on the results, and the gather the log-likelihoods into two columns
  # named "inferred_collection" and "log_likelihood"
  result <- reference %>%
    dplyr::select(indiv, collection, repunit) %>%
    dplyr::bind_cols(., tibble::as_tibble(logl)) %>%
    tidyr::gather(data = ., key = "inferred_collection", value = "log_likelihood", -indiv, -collection, -repunit) %>%
    dplyr::mutate(indiv = factor(indiv, levels = unique(indiv))) %>% # this lets us keep indivs in input order
    dplyr::arrange(indiv, desc(log_likelihood)) %>%
    dplyr::mutate(indiv = as.character(indiv))  # when done, coerce indiv back to character

  # we need to do something similar with the sums-of-squares of the logls, then
  # we will left_join it
  ssq_tibble <- reference %>%
    dplyr::select(indiv, collection, repunit) %>%
    dplyr::bind_cols(., tibble::as_tibble(logl_ssq)) %>%
    tidyr::gather(data = ., key = "inferred_collection", value = "ssq_logl", -indiv, -collection, -repunit)

  # here we join that on
  result <- dplyr::left_join(result, ssq_tibble,
                      by = c("indiv", "collection", "repunit", "inferred_collection"))

  # and finally, we use a join to put a column on there for "inferred_repunit".
  # this ugly thing just gets a tibble that associates repunits with collections
  repu_assoc <- result %>%
    dplyr::count(collection, repunit) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup() %>%
    dplyr::rename(inferred_collection = collection,
                  inferred_repunit = repunit)

  # and this joins the inferred_repunits column on there and then
  # orders the columns in a good way, and finally adds a column of
  # scaled likelihoods for each individual, then ungroups and left_joins
  # to the number of loci.
  result %>%
    dplyr::left_join(., repu_assoc, by = "inferred_collection") %>%
    dplyr::select(indiv:inferred_collection, inferred_repunit, log_likelihood, ssq_logl) %>%
    dplyr::group_by(indiv) %>%
    dplyr::mutate(scaled_likelihood = exp(log_likelihood) / sum(exp(log_likelihood))) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., count_missing_data(reference, gen_start_col), by = "indiv")


}

