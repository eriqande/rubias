#include <Rcpp.h>
using namespace Rcpp;

double rgammadouble(int a, double b, double c);
NumericVector dirch_from_allocations(IntegerVector C, NumericVector lambda);
NumericVector dirch_from_counts(IntegerVector C, NumericVector lambda);

int randint(int n);
IntegerVector tabulate_allocations(IntegerVector C, int n);
IntegerVector turn_non_catch_to_minus_one(IntegerVector Z, NumericVector P, int &NC);
