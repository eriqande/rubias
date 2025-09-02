#include <Rcpp.h>
using namespace Rcpp;


//' Sample 1 observation from cell probabilities that are columns of a matrix
//'
//' Takes a matrix in which columns sum to one. For each column, performs a
//' single multinomial draw from the rows, weighted by their values in that column
//'
//' @keywords internal
//'
//' @param M a matrix whose columns are reals summing to one
//'
//' @return a vector length = \code{ncol(M)} of indices, with each element being
//' the row that was chosen in that column's sampling
//' @export
// [[Rcpp::export]]
IntegerVector samp_from_mat(NumericMatrix M) {
  int C = M.ncol();
  int R = M.nrow();
  int r,c, res;
  double cumul, rando;
  IntegerVector out(C);
  NumericVector rando_vec(C);

  rando_vec = runif(C);

  for(c = 0; c < C; c++) {
    cumul = 0.0;
    rando = rando_vec[c];
    for(r = 0; r < R; r++) {
      cumul += M(r, c);
      res = r + 1L;
      if(cumul >= rando) {
        break;
      }
    }
    out[c] = res;
  }
  return(out);
}


//' Simulate a single multinomial vector from within Rcpp
//'
//' From: https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/
//'
//' @keywords internal
//' @param size the number of trials
//' @param probs  the cell probabilities
//' @param N the number of cells
//' @return An IntegerVector of length N
//' @export
// [[Rcpp::export]]
IntegerVector rmultinom_1(unsigned int &size, NumericVector &probs, unsigned int &N) {
  IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}
