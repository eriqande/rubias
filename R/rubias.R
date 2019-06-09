#' rubias: Bayesian inference from the conditional genetic stock identification model
#'
#' Read the "rubias-overview" vignette for information on data input formats
#' and how to use the package
#'
#'
#' @section the \code{rubias} main high-level functions:
#'
#' The following functions are wrappers, designed for user-friendly input and useful output:
#'
#' \code{\link{infer_mixture}} is used to perform genetic stock identification.
#' Options include standard MCMC and the parametric bootstrap bias correction.
#'
#' \code{\link{self_assign}} does simple self-assignment of individuals in a reference data set
#' to the various collections in the reference data set.
#'
#' \code{\link{assess_reference_loo}} does leave-one-out based simulations to predict how
#' accurately GSI can be done.
#'
#' \code{\link{assess_reference_mc}} uses Monte-Carlo cross-validation based simulations
#' to predict how accurately GSI can be done.
#'
#' \code{\link{assess_pb_bias_correction}} attempts to demonstrate how much (or little)
#' improvement can be expected from the parametric bootstrap correction given a particular
#' reference data set.
#'
#'
#' @section genetic data format:
#'
#' See the vignette.
#'
#' @section example data:
#'
#' \code{\link{alewife}}, \code{\link{blueback}}, and \code{\link{chinook}} are
#' genetic data sets that are useful for playing around with rubias and testing it
#' out.
#'
#' @docType package
#' @name rubias
#' @importFrom stats rbeta rmultinom var
#' @importFrom utils write.table
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib rubias
NULL


# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "allele",
      "bh_rho",
      "bias",
      "coll_int",
      "collection",
      "collection_1",
      "collection_2",
      "collection_scenario",
      "dev",
      "expected_mean",
      "expected_var",
      "indiv",
      "indiv_1",
      "indiv_2",
      "indx1",
      "indx2",
      "inferred_collection",
      "inferred_repunit",
      "iter",
      "locus",
      "log_likelihood",
      "mc_pi",
      "mean_bias",
      "mean_prop_bias",
      "method",
      "mle",
      "mse",
      "n",
      "normo_logl",
      "num_match",
      "num_non_miss",
      "omega",
      "pb_rho",
      "pofz",
      "post_mean",
      "prop_bias",
      "repunit",
      "repunit_scenario",
      "rho",
      "rho_est",
      "rho_mcmc",
      "rho_pb",
      "sample_type",
      "scaled_likelihood",
      "simulated_collection",
      "simulated_repunit",
      "times_seen",
      "true_rho"
    )
  )
}


