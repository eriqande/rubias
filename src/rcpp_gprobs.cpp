#include <Rcpp.h>
#include "macros.h"
using namespace Rcpp;

//' Draw all allele counts for a population from a vector of allele counts
//'
//' (This is a test for our allele count vector indexing macro)
//' @param c the collection index desired
//' @param l the locus desired
//' @param ac_v the vector of allele counts
//' @param L the total number of loci
//' @param C the total number of collections
//' @param A an integer vector of alleles in each population
//' @param CA an integer vector of the cumulative number of alleles at all previous loci in the vector
//'
//' @export
// [[Rcpp::export]]
IntegerVector coll_ac(int l, int c, IntegerVector ac_v, int L, int C, IntegerVector A, IntegerVector CA) {
  IntegerVector out(A[l]);
  int a;
  for(a = 0; a < A[l]; a++) {
    out[a] = ac_v[D_dx(l, c, a, L, C, A, CA)];
  }
  return out;
}

//' Fetch individual genotypes from a locus and population
//'
//' (This is a test case for our individual indexing macro)
//'
//' @param i the index of the desired individual to sample
//' @param l the index of the desired locus
//' @param I_v a vector of individual genotypes, created by unlisting allelic_list output
//' @param P Ploidy
//' @param I Total number of individuals
//'
//' @export
// [[Rcpp::export]]
IntegerVector genotype_i(int l, int i, IntegerVector I_v, int P, int I) {
  IntegerVector out(P);
  int g;
  for(g = 0; g < P; g++) {
    out[g] = I_v[I_dx(l, i, g, P, I)];
  }
  return out;
}

// test/example for the Dirichlet indexing macro
/*** R

example(a_freq_list)
ale_ac_v <- unlist(ale_ac)
ale_ac$Ap058[ ,"ROA"]
A <- sapply(ale_ac, nrow)
CA <- cumsum(c(0,A))
coll_ac(9, 17, ale_ac_v, length(A), length(ale_ac$Ap058), A, CA)
*/
// test/example for the individual genotype indexing macro
/*** R

example(allelic_list)
ale_alle_list <- lapply(ale_alle_list, function(l) {
  t(cbind(l$a, l$b))
})
ale_alle_v <- unlist(ale_alle_list)
ale_alle_list$Aa091[,3]
genotype_i(4,2, ale_alle_v, 2, 1070)
*/


//' Calculate list of genotype likelihoods
//'
//' @param par_list genetic data converted to the param_list format by tcf2param_list
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

  for(i = 0; i < N; i++) {
    for(c = 0; c < C; c++) {
      sum = 0.0;
      LOO = c == (coll[i] - 1);
      for(l = 0; l < L; l++) {
        GPROB_DIP(i, l, c, gp);
        sum += log(gp);
      }
      out(c, i) = sum;
    }
  }
  return(out);
}

//' Calculate list of genotype likelihoods
//'
//' @param par_list genetic data converted to the param_list format by tcf2param_list
//'
//' @examples
//' example(tcf2param_list)
//' ale_glL <- geno_logL(ale_par_list)
//' @export
// [[Rcpp::export]]
NumericMatrix geno_logL_RU(List par_list) {
  int i, c, l;
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  int N = as<int>(par_list["N"]);
  int C = as<int>(par_list["K"]);
  int L = as<int>(par_list["L"]);
  int LOO;
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  IntegerVector ru = as<IntegerVector>(par_list["ru"]);
  NumericVector DP = as<NumericVector>(par_list["DP_R"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP_R"]);
  double sum, gp;
  NumericMatrix out(C, N);

  for(i = 0; i < N; i++) {
    for(c = 0; c < C; c++) {
      sum = 0.0;
      LOO = c == (ru[i] - 1);
      for(l = 0; l < L; l++) {
        GPROB_DIP(i, l, c, gp);
        sum += log(gp);
      }
      out(c, i) = sum;
    }
  }
  return(out);
}
