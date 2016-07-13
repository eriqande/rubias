#include <Rcpp.h>
using namespace Rcpp;


double rgammadouble(int a, double b, double c)
{   Rcpp::NumericVector x = rgamma(a,b,1/c);
  return x(0);
}

// Given a maximum integer n, simulates an integer from 1:n
int randint(int n) {
  if(n == 0) {
    return(0);
  }
  int ind;
  double rando;
  rando = runif(1)[0];
  ind = ceil(rando * n);
  return(ind);
}

// Some simple utilities for Rcpp

//' given a vector of different categories in 1...n and a prior simulate a dirichlet r.v.
//'
//' The categories are labeled in C from 1 up to n.  n is the length of pi which is vector of priors.
//' Note that all elements of pi must be strictly greater than 0.
//' @param C  a vector giving different categories of individual (not counts of categories---must be tabulated)
//' @param pi  priors for the categories
//' @export
// [[Rcpp::export]]
NumericVector dirch_from_allocations(IntegerVector C, NumericVector lambda) {
  int i;
  int n = lambda.size();
  int N = C.size();
  NumericVector out = clone(lambda);

  for(i = 0; i < N; i++) {
    out[C[i] - 1] += 1.0;
  }
  for(i = 0; i < n; i++) {
    out[i] = rgammadouble(1L, out[i], 1.0);
  }

  return(out / sum(out));
}

// Some simple utilities for Rcpp

//' given a vector of different categories in 1...n and a prior simulate a dirichlet r.v.
//'
//' The categories are labeled in C from 1 up to n.  n is the length of pi which is vector of priors.
//' Note that all elements of pi must be strictly greater than 0.
//' @param C  a vector giving counts of categories
//' @param pi  priors for the categories
//' @export
// [[Rcpp::export]]
NumericVector dirch_from_counts(IntegerVector C, NumericVector lambda) {
  int i;
  int n = lambda.size();
  NumericVector out = clone(lambda);

  for(i = 0; i < n; i++) {
    out[i] += C[i];
  }
  for(i = 0; i < n; i++) {
    out[i] = rgammadouble(1L, out[i], 1.0);
  }

  return(out / sum(out));
}
