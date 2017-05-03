#include <Rcpp.h>
#include "macros.h"
using namespace Rcpp;

//' Calculate a matrix of genotype log-likelihoods for a genetic dataset
//'
//' Takes a list of key parameters from a genetic dataset, and calculates
//' the log-likelihood of each individual's genotype, given the allele counts
//' in each collection
//'
//' Leave-One-Out cross-validation is used to avoid bias in log-likelihood for an
//' individual's known collection of origin
//'
//' @keywords internal
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//'
//' @return \code{geno_logL} returns a matrix with C rows and I columns. Each column
//' represents an individual, and each row the log-likelihood of that individual's
//' genotype, given the allele counts in that collection
//'
//' @examples
//' example(tcf2param_list)
//' ale_glL <- geno_logL(ale_par_list)
//' @export
// [[Rcpp::export]]
NumericMatrix geno_logL(List par_list) {
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
  double sum, gp;
  NumericMatrix out(C, N);

  for(i = 0; i < N; i++) { // cycle over individuals
    for(c = 0; c < C; c++) { // cycle over collections
      sum = 0.0;
      LOO = c == (coll[i] - 1);
      for(l = 0; l < L; l++) {  // cycle over loci
        GPROB_DIP(i, l, c, gp);
        sum += log(gp);
      }
      out(c, i) = sum;
    }
  }
  return(out);
}
