#include <Rcpp.h>
using namespace Rcpp;
#include "utilities.h"
#include "rcpp_sampling.h"


//' EM algorithm from the simplest GSI model for pi and the individual posterior probabilities
//'
//' Using a matrix of scaled likelihoods, this function does an EM algorithm to climb the
//' likelihood surface for pi, and computes the plug-in estimate of the posteriors
//' for all the individuals.  It returns the output in a list.
//' @keywords internal
//' @param SL  a matrix of the scaled likelihoods.  This is should have values for each individual in a column
//' (going down in the rows are values for different collections).
//' @param Pi_init  Starting value for the pi (collection mixture proportion) vector.
//' @param max_iterations the maximum total number of reps iterations to do.
//' @param tolerance the EM-algorithm will be considered converged when the sum over the elements of pi of the absolute value
//' of the difference between the previous and the current estimate is less than tolerance.
//' @param return_progression  If true, then the pi_trace component of the output shows the value of pi visited en route to the end.
//'
//' @return \code{gsi_em_1} returns a final Maximum-Likelihood estimate for pi and PofZ,
//' as well as the number of iterations needed to reach convergence ("iterations_performed"),
//' and traces of the pi values and change in pi in each iteration
//'
//' @examples
//' params <- tcf2param_list(alewife, 15)
//' logl <- geno_logL(params)
//' SL <- apply(exp(logl), 2, function(x) x/sum(x))
//' test_em <- gsi_em_1(SL, rep(1/params$C, params$C), max_iterations = 10^6, tolerance = 10^-7, return_progression = TRUE)
//'
//' @export
// [[Rcpp::export]]
List gsi_em_1(NumericMatrix SL, NumericVector Pi_init, int max_iterations, double tolerance, bool return_progression) {
  List pi_list;
  List ret;
  NumericVector pi = clone(Pi_init);
  NumericVector new_pi = clone(Pi_init);
  NumericVector diff_vec(max_iterations, NA_REAL);
  NumericMatrix posts = clone(SL);

  int R = SL.nrow();
  int C = SL.ncol();
  int i, r, c;
  double sum, tmp, diff;

  // store pi value
  if( return_progression == true ) {
    pi_list.push_back(pi);
  }

  for(i = 0; i < max_iterations; i++) {

    // initialize new_pi to 0.0's to prepare for accumulating a sum
    for(r = 0; r < R; r++) { new_pi[r] = 0.0; }

    // normalize the scaled likelihoods into posteriors
    for(c = 0; c < C; c++) {
      sum = 0.0;
      for(r = 0; r < R; r++) {
        tmp = SL(r, c) * pi[r];
        posts(r, c) = tmp;
        sum += tmp;
      }
      for(r = 0; r < R; r++) {
        posts(r, c) /= sum;
        new_pi[r] += posts(r, c) / (double)C;
      }
    }

    // divide new_pi by C to get the means
//    new_pi = new_pi / C;

    // compute the sum of absolute differences then store it
    diff = 0.0;
    for(r = 0; r < R; r++) {  // I was unable to get this to compile when using equivalent Rcpp sugar functions.
      diff += (new_pi[r] - pi[r]) * (1 - 2 * ((new_pi[r] - pi[r]) < 0.0));
      pi[r] = new_pi[r];

    }
    diff_vec[i] = diff;

    // store pi value if requested
    if( return_progression == true ) {
      pi_list.push_back(pi);
    }

    if(diff < tolerance) {
       break;  // break out if there is little change from the last value.
    }
  }

  // note that we use the posteriors from the next to last iteration of pi, but that should be fine...
  ret = List::create(pi, posts, i, diff_vec, pi_list);
  ret.names() = CharacterVector::create("pi", "PofZ", "iterations_performed", "diffs_vector", "pi_trace");
  return(ret);
}
