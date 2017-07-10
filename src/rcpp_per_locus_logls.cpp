#include <Rcpp.h>
#include "macros.h"
using namespace Rcpp;

//' Return a matrix of locus-specific self-assignment logls
//'
//' Takes a list of key parameters from a genetic dataset, and calculates
//' the log-likelihood of each individual's single-locus genotype, given the allele counts
//' in the individual's collection.
//'
//' This uses Leave-One-Out cross-validation is used to avoid bias in log-likelihood for an
//' individual's known collection of origin
//'
//' @keywords internal
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//'
//' @return \code{rcpp_per_locus_logls} returns a matrix with I rows and L columns. Each row
//' represents an individual, and each column a locus. Note that missing data at a locus
//' returns a zero.  That should be changed to NA later.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix rcpp_per_locus_logls(List par_list) {
  int i, c, l;
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  int N = as<int>(par_list["N"]);
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);
  int LOO;
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  IntegerVector coll = as<IntegerVector>(par_list["coll"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  double gp;
  NumericMatrix out(N, L);

  for(i = 0; i < N; i++) { // cycle over individuals
    LOO = 1;
    c = (coll[i] - 1);  // only compute logl for the population the individual is from
    for(l = 0; l < L; l++) {  // cycle over loci
      GPROB_DIP(i, l, c, gp);
      out(i, l) = log(gp);
    }
  }
  return(out);
}
