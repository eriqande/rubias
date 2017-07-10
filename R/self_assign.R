#' Do leave-one-out self-assignment of individuals in a reference baseline
#'
#' Returns a tidy data frame
#' @inheritParams assess_reference_loo
#' @param preCompiledParams Users should never use this option.  It is here only so that
#' this function can be called on a precompiled set of parameters with infer_mixture.  Don't
#' use this, unless you are one of the package developers...
#' @return a tibble ...
#' @export
#' @examples
#' ale_sa <- self_assign(alewife, 17)
self_assign <- function(reference, gen_start_col, preCompiledParams = NULL) {

  # if not supplying it with preCompiledParams, check reference and get the params
  if (is.null(preCompiledParams)) {
    # make sure the reference file is OK
    check_refmix(reference, gen_start_col, "reference")

    # get the necessary parameters from the reference data
    params <- tcf2param_list(reference, gen_start_col, summ = T)
  } else {
    # otherwise, assume the reference was checked elsewhere (i.e. in infer_mixture())
    # and just set params to preCompiledParams
    params <- preCompiledParams
  }

  # get the log-likelihoods
  logl <- t(geno_logL(par_list = params))
  # and get the sum-of-squares over loci of the logls (for z-score calculation)
  logl_ssq <- t(geno_logL_ssq(par_list = params))

  # put the collection names at the top of them. To do this, we put RU_vec into sorted
  # order and then grab the names off it
  colnames(logl) <- names(sort(params$RU_vec))
  colnames(logl_ssq) <- names(sort(params$RU_vec))

  # for locus-specific logls calculation, get the locus-specific means and vars for each collection
  locus_means_and_vars <- per_locus_means_and_vars(params)
  # then compute exptected means and vars for each individual on the basis of patterns of missing data
  mv_sums <- rcpp_indiv_specific_logl_means_and_vars(params, locus_means_and_vars)


  # then get a tibble with indiv, inferred_collection, log_likelihood, and the expected mean and var
  logl_tibble <- tibble::tibble(
    indiv = rep(params$indiv_names, each = params$C),
    inferred_collection = rep(params$collection_names, params$N),
    log_likelihood = as.vector(t(logl)),
    expected_mean = as.vector(mv_sums$mean),
    expected_var = as.vector(mv_sums$var)
  )


  # then make a tibble of the meta data (indiv, collection, repuunit) from
  # "reference" back on the results, join on the logls and expected means and variances of the logls,
  # sort by logl with indivs in input order, and then compute the z_score
  result <- reference %>%
    dplyr::select(indiv, collection, repunit) %>%
    dplyr::left_join(logl_tibble, by = "indiv") %>%
    dplyr::mutate(indiv = factor(indiv, levels = unique(indiv))) %>% # this lets us keep indivs in input order
    dplyr::arrange(indiv, desc(log_likelihood)) %>%
    dplyr::mutate(indiv = as.character(indiv))  %>% # when done, coerce indiv back to character
    dplyr::mutate(z_score = (log_likelihood - expected_mean) / sqrt(expected_var)) %>%  # compute z_score
    dplyr::select(-expected_mean, -expected_var)  # remove the columns used for the z_score calculation



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
  ret <- result %>%
    dplyr::left_join(., repu_assoc, by = "inferred_collection") %>%
    dplyr::group_by(indiv) %>%
    dplyr::mutate(normo_logl = log_likelihood - mean(log_likelihood)) %>%  # This is done to prevent underflow
    dplyr::mutate(scaled_likelihood = exp(normo_logl) / sum(exp(normo_logl))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-normo_logl) %>%
    dplyr::left_join(., count_missing_data(reference, gen_start_col), by = "indiv") %>%
    dplyr::select(indiv:inferred_collection, inferred_repunit, scaled_likelihood, dplyr::everything())


  ret
}

