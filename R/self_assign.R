#' Do leave-one-out self-assignment of individuals in a reference baseline
#'
#' Returns a tidy data frame
#' @inheritParams assess_reference_loo
#' @return a tibble ...
#' @export
self_assign <- function(reference, gen_start_col) {
  # get the necessary parameters from the reference data
  params <- tcf2param_list(reference, gen_start_col, summ = T)

  # get the log-likelihoods
  logl <- t(geno_logL(par_list = params))

  # put the collection names at the top of them. To do this, we put RU_vec into sorted
  # order and then grab the names off it
  colnames(logl) <- names(sort(params$RU_vec))

  # then make a tibble of it and put the meta data (indiv, collection, repuunit) from
  # "reference" back on the results, and the gather the log-likelihoods into two columns
  # named "inferred_collection" and "log_likelihood"
  result <- reference %>%
    dplyr::select(indiv, collection, repunit) %>%
    dplyr::bind_cols(., tibble::as_tibble(logl)) %>%
    tidyr::gather(data = ., key = "inferred_collection", value = "log_likelihood", -indiv, -collection, -repunit) %>%
    dplyr::arrange(indiv, desc(log_likelihood))

  # then, if collection is a factor, we make inferred_collection a factor with the same levels
  if (is.factor(result$collection)) {
    result$inferred_collection <- factor(result$inferred_collection,
                                         levels = levels(result$collection))
  }

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
  # scaled likelihoods for each individual, then ungroups and returns
  # the result
  result %>%
    dplyr::left_join(., repu_assoc, by = "inferred_collection") %>%
    dplyr::select(indiv:inferred_collection, inferred_repunit, log_likelihood) %>%
    dplyr::group_by(indiv) %>%
    dplyr::mutate(scaled_likelihood = exp(log_likelihood) / sum(exp(log_likelihood))) %>%
    dplyr::ungroup()
}

