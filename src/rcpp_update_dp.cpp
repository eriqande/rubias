#include <Rcpp.h>
#include "macros.h"
using namespace Rcpp;

//' Update a matrix of Dirichlet parameters given a set of collection allocations
//'
//' Takes a list of key parameters from a genetic dataset, as well as a vector
//' of collections to which individuals have been allocated, and updates the baseline
//' allele frequencies as if the individuals were known to originate from said collections
//'
//' For each individual and locus, this identifies the alleles held by that individual,
//' then adds 1 to these alleles' counts in the Dirichlet parameter vector
//'
//' @keywords internal
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//' @param allocations an N-long vector of collection ID numbers from \code{samp_from_post}
//'
//' @return \code{geno_logL} returns a matrix with C rows and I columns. Each column
//' represents an individual, and each row the log-likelihood of that individual's
//' genotype, given the allele counts in that collection
//' @export
// [[Rcpp::export]]
NumericVector update_dp(List par_list, IntegerVector allocations) {
  int i, c, l;
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  int N = as<int>(par_list["N"]);
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  IntegerVector PLOID = as<IntegerVector>(par_list["ploidies"]);
  NumericVector DP_temp = clone(DP);

  for(i = 0; i < N; i++) {
    for(l = 0; l < L; l++) {
      int a1 = I[I_dx(l, i, 0, 2, N)] - 1;
      int a2 = I[I_dx(l, i, 1, 2, N)] - 1;
      c = allocations[i] - 1;

      if(PLOID[l] == 1) {
        if(a1 >= 0) {
          DP_temp[D_dx(l, c, a1, L, C, A, CA)] += 1;
        }
      } else {
        if(a1 >= 0 && a2 >= 0) {
          DP_temp[D_dx(l, c, a1, L, C, A, CA)] += 1;
          DP_temp[D_dx(l, c, a2, L, C, A, CA)] += 1;
        }
      }
    }
  }
  return(DP_temp);
}
