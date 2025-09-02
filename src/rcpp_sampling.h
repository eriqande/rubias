#include <Rcpp.h>
using namespace Rcpp;

IntegerVector samp_from_mat(NumericMatrix M) ;
IntegerVector rmultinom_1(unsigned int &size, NumericVector &probs, unsigned int &N);
