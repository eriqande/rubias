#include <Rcpp.h>
using namespace Rcpp;


//' Sample 1 observation from cell probabilities that are columns of a matrix
//'
//' blah blah...
//' @param M a matrix whose columns are reals summing to one
//'
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
