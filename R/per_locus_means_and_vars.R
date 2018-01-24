

#' Compute the mean and variance of the single-locus genotype likelihoods for each collection
#'
#' This assumes that you have compiled params for a reference data set and then it just
#' calls rcpp_per_locus and then summarizes the results.
#' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}. This should be
#' include genotypes only for the reference individuals.
#' @return Returns a list with two components, mean and var, each one a matrix that has
#' C (number of collections) rows and L (number of loci) columns, giving the mean (or variance)
#' of the genotype likelihoods in the individuals in that collection at that locus.
#' @keywords internal
#' @export
per_locus_means_and_vars <- function(par_list) {

  gl <- rcpp_per_locus_logls(par_list)  # get the matrix of logls

  # turn the zeroes into NAs.  Note that with the allele frequency priors we use, we
  # never expect to have a genotype probability of 1.  Hence we know that log-geno-probs
  # of zero must be missing loci.
  gl[dplyr::near(gl, 0.0)] <- NA

  gld <- as.data.frame(gl)

  # Now, for each collection we want to compute means and variances.  We will do this by splitting
  # it and then lapplying and binding rows at the end
  gl_list <- split(gld, par_list$coll)

  # then compute the means
  gl_means <- lapply(gl_list, function(x) colMeans(as.matrix(x), na.rm = TRUE)) %>%
    do.call(args = ., what = rbind)

  # and compute the variances
  gl_vars <- lapply(gl_list, function(x) {
    apply(x, 2, var, na.rm = TRUE)
  }) %>%
    do.call(args = ., what = rbind)

  # now, there is an awkward thing where every individual in a collection might be missing
  # data at a given locus.  These elements are NaNs in gl_means.  We will turn their menas
  # and also their variances to zero.
  gl_means[is.na(gl_means)] <- 0.0
  gl_vars[is.na(gl_vars)] <- 0.0

  list(mean = gl_means, var = gl_vars)
}
