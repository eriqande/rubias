#include <Rcpp.h>
#include "macros.h"
using namespace Rcpp;

//' From the pattern of missing data at each individual, compute the expected mean and variance of the logl
//'
//' This takes a param_list so that it has access to individual genotypes (and hence can cycle through them
//' and know which are missing and which are not.)  It also takes a matrix of per-locus logl means
//' and variances like what is computed by \code{\link{per_locus_means_and_vars}}.
//'
//' This function doesn't do any checking to assure that the par_list and the per-locus logl means
//' matrix are made for one another.  (i.e. use the same collections in the same order.)
//'
//' @keywords internal
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//' @param MV a list of two matrices, one of means and the other of variances, which are C x L
//' matrices.  This is basically the list that is returned by \code{\link{per_locus_means_and_vars}}.
//'
//' @return a matrix with C rows and I columns. Each row
//' represents a collection, and each column an individual.
//'
//' @export
// [[Rcpp::export]]
List rcpp_indiv_specific_logl_means_and_vars(List par_list, List MV) {
  int i, c, l;
  double meansum, varsum, miss_flag;
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  int N = as<int>(par_list["N"]);
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);

  NumericMatrix means = as<NumericMatrix>(MV["mean"]);
  NumericMatrix vars = as<NumericMatrix>(MV["var"]);

  NumericMatrix outmean(C, N);
  NumericMatrix outvar(C, N);

  for(i = 0; i < N; i++) { // cycle over individuals
    for(c = 0; c < C; c++) {  // cycle over collections
      meansum = 0.0; // initialize to accumulate a sum over non-missing loci in the individual
      varsum = 0.0;

      for(l = 0; l < L; l++) {  // cycle over loci
        ZERO_IF_LOCUS_MISSING(i, l, c, miss_flag);
        meansum += miss_flag * means(c, l);
        varsum += miss_flag * vars(c, l);
      }

      outmean(c, i) = meansum;
      outvar(c, i) = varsum;
    }
  }
  return List::create(Named("mean") = outmean,
                      Named("var") = outvar);
}
