#include <Rcpp.h>
using namespace Rcpp;
#include "utilities.h"
#include "rcpp_sampling.h"

//' MCMC from the simplest GSI model for pi and the individual posterior probabilities
//'
//' Using a matrix of scaled likelihoods, this function samples values of pi and the posteriors
//' for all the individuals.  It returns the output in a list.
//' @keywords internal
//' @param SL  matrix of the scaled likelihoods.  This is should have values for each individual in a column
//' (going down in the rows are values for different populations).
//' @param Pi_init  Starting value for the pi (collection mixture proportion) vector.
//' @param lambda the prior to be added to the collection allocations, in order to
//' generate pseudo-count Dirichlet parameters for the simulation of a new pi vector
//' @param reps total number of reps (sweeps) to do.
//' @param burn_in how many reps to discard in the beginning when doing the mean calculation. They will still be
//' returned in the traces if desired
//' @param sample_int_Pi the number of reps between samples being taken for Pi traces.  If 0 no trace samples are taken
//' @param sample_int_PofZ the number of reps between samples being taken for the traces of posterior of each individual's origin. If 0
//' no trace samples are taken.
//'
//' @return \code{gsi_mcmc_1} returns a list of three. \code{$mean} lists the posterior
//' means for collection proportions \code{pi} and for the individual posterior
//' probabilities of assignment \code{PofZ}. \code{$sd} returns the posterior standard
//' deviations for the same values.
//'
//' If the corresponding \code{sample_int} variables are not 0, \code{$trace} contains
//' samples taken from the Markov chain at intervals of \code{sample_int_}(variable) steps.
//'
//' @examples
//' # this demonstrates it with scaled likelihoods computed from
//' # assignment of the reference samples
//' params <- tcf2param_list(alewife, 17)
//' logl <- geno_logL(params)
//' SL <- apply(exp(logl), 2, function(x) x/sum(x))
//' lambda <- rep(1/params$C, params$C)
//' mcmc <- gsi_mcmc_1(SL, lambda, lambda, 200, 50, 5, 5)
//' @export
// [[Rcpp::export]]
List gsi_mcmc_1(NumericMatrix SL, NumericVector Pi_init, NumericVector lambda, int reps, int burn_in, int sample_int_Pi, int sample_int_PofZ) {
  List pi_list;
  List PofZ_list;
  List trace, mean, sd, ret;
  NumericVector pi = clone(Pi_init);
  NumericVector pi_sums(Pi_init.size());
  NumericVector pi_sums_sq(Pi_init.size());
  NumericMatrix posts = clone(SL);
  NumericMatrix post_sums(SL.nrow(), SL.ncol());
  NumericMatrix post_sums_sq(SL.nrow(), SL.ncol());
  NumericMatrix sd_ret(SL.nrow(), SL.ncol());

  int R = SL.nrow();
  int C = SL.ncol();
  int i, r, c;
  double sum, tmp;
  int num_samp = reps - burn_in;

  if(num_samp <= 1) stop("reps - burn_in <= 1");

  for(i = 0; i < reps; i++) {
    // store pi value
    if( (sample_int_Pi > 0) && (i % sample_int_Pi == 0) ) {
      pi_list.push_back(pi);
    }
    if(i >= burn_in) {
      pi_sums += pi;
      pi_sums_sq += pi * pi;
    }
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
        if(i >= burn_in) {
          post_sums(r, c) += posts(r, c);
          post_sums_sq(r, c) += posts(r, c) * posts(r, c);
        }
      }
    }

    // store PofZ values
    if( (sample_int_PofZ > 0) && (i % sample_int_PofZ == 0) ) {
      PofZ_list.push_back(posts);
    }

    // allocate individuals to populations and simulate a new pi
    pi = dirch_from_allocations(samp_from_mat(posts), lambda);

  }

  // put the traces in there if there are any
  trace = List::create(pi_list, PofZ_list);
  trace.names() = CharacterVector::create("pi", "PofZ");

  // put the means and standard devs and traces in the return variable'
  post_sums = post_sums / num_samp;
  pi_sums  = pi_sums / num_samp;

  mean.push_back(pi_sums);
  sd.push_back(sqrt((pi_sums_sq - (num_samp * pi_sums * pi_sums)) / (num_samp - 1.0)));

  mean.push_back(post_sums);
  // had to write this without Rcpp sugar; seems to self-destruct trying to run on a matrix
  for(c = 0; c < C; c++){
    for(r = 0; r < R; r++){
      sd_ret(r,c) = sqrt((post_sums_sq(r,c) - (num_samp * post_sums(r,c) * post_sums(r,c))) / (num_samp - 1.0));
    }
  }
  sd.push_back(sd_ret);

  mean.names() = CharacterVector::create("pi", "PofZ");
  sd.names() = CharacterVector::create("pi", "PofZ");

  ret = List::create(mean, sd, trace);
  ret.names() = CharacterVector::create("mean", "sd", "trace");
  return(ret);
}
