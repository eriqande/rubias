#include <Rcpp.h>
#include "utilities.h"
#include "macros.h"
#include "rcpp_sampling.h"
#include "rcpp_update_dp.h"
using namespace Rcpp;

//' MCMC from the fully Bayesian GSI model for pi and the individual posterior probabilities
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
//' @return \code{gsi_mcmc_fb} returns a list of three. \code{$mean} lists the posterior
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
//'
//' # we have to get the ploidies to pass to tcf2param_list
//' locnames <- names(alewife)[-(1:16)][c(TRUE, FALSE)]
//' ploidies <- rep(2, length(locnames))
//' names(ploidies) <- locnames
//'
//' params <- tcf2param_list(alewife, 17, ploidies = ploidies)
//' logl <- geno_logL(params)
//' SL <- apply(exp(logl), 2, function(x) x/sum(x))
//' lambda <- rep(1/params$C, params$C)
//' mcmc <- gsi_mcmc_1(SL, lambda, lambda, 200, 50, 5, 5)
//' @export
// [[Rcpp::export]]
List gsi_mcmc_fb(List par_list, NumericVector Pi_init, NumericVector lambda,
                 int reps, int burn_in, int sample_int_Pi, int sample_int_PofZ) {

  // Code from geno_logl, creating C X N loglikelihood matrix
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
  IntegerVector PLOID = as<IntegerVector>(par_list["ploidies"]);
  double sum, gp;
  NumericVector colsums(N);
  NumericVector colmeans(N);
  NumericMatrix logl(C, N);
  NumericMatrix sweep_logl(C, N);
  NumericMatrix SL(C, N);

  // code from gsi_mcmc_1
  int r;
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
  IntegerVector allocs(SL.ncol());
  NumericVector DP_temp(DP.size());

  double tmp;
  int num_samp = reps - burn_in;

  for(r = 0; r < reps; r++) {
    // store pi value
    if( (sample_int_Pi > 0) && (r % sample_int_Pi == 0) ) {
      pi_list.push_back(pi);
    }
    if(r >= burn_in) {
      pi_sums += pi;
      pi_sums_sq += pi * pi;
    }

    for(i = 0; i < N; i++) { // cycle over individuals
    for(c = 0; c < C; c++) { // cycle over collections
      sum = 0.0;
      LOO = c == (coll[i] - 1);
      for(l = 0; l < L; l++) {  // cycle over loci
        GPROB(i, l, c, gp);
        sum += log(gp);
      }
      logl(c, i) = sum;
      colsums(i) += sum; // sum across collections for column mean calculation
    }
  }

  // take column means, then sweep out from each logl
  for(i = 0; i < N; i++) {
    colmeans(i) = colsums(i)/C;
    for(c = 0; c < C; c++) {
      sweep_logl(c, i) = logl(c, i) - colmeans(i);
    }
  }

  // convert to scaled logl matrix SL
  for(i = 0; i < N; i++) {
    sum = 0.0;
    for(c = 0; c < C; c++) {
      sum += exp(sweep_logl(c, i));
    }
    for(c = 0; c < C; c++) {
      SL(c, i) = exp(sweep_logl(c, i)) / sum;
    }
  }

  // code from gsi_mcmc_1



  if(num_samp <= 1) stop("reps - burn_in <= 1");


    // normalize the scaled likelihoods into posteriors
    for(i = 0; i < N; i++) {
      sum = 0.0;
      for(c = 0; c < C; c++) {
        tmp = SL(c, i) * pi[c];
        posts(c, i) = tmp;
        sum += tmp;
      }
      for(c = 0; c < C; c++) {
        posts(c, i) /= sum;
        if(r >= burn_in) {
          post_sums(c, i) += posts(c, i);
          post_sums_sq(c, i) += posts(c, i) * posts(c, i);
        }
      }
    }

    // store PofZ values
    if( (sample_int_PofZ > 0) && (r % sample_int_PofZ == 0) ) {
      PofZ_list.push_back(posts);
    }

    // allocate individuals to populations and simulate a new pi
    allocs = samp_from_mat(posts);
    pi = dirch_from_allocations(allocs, lambda);
    // compute a new Dirichlet Parameter Vector based on the allocations
    DP = update_dp(par_list, allocs);

    for(i = 0; i < N; i++) {
      c = allocs[i] - 1;
      for(l = 0; l < L; l++) {
        if(PLOID[l] == 1) {
          sum_DP[SD_dx(l, c, C)] += 1;
        } else {
          sum_DP[SD_dx(l, c, C)] += 2;
        }
      }
    }
  }

  // put the traces in there if there are any
  trace = List::create(pi_list, PofZ_list);
  trace.names() = CharacterVector::create("pi", "PofZ");

  // put the means and standard devs and traces in the return variable
  post_sums = post_sums / num_samp;
  pi_sums  = pi_sums / num_samp;

  mean.push_back(pi_sums);
  sd.push_back(sqrt((pi_sums_sq - (num_samp * pi_sums * pi_sums)) / (num_samp - 1.0)));

  mean.push_back(post_sums);
  // had to write this without Rcpp sugar; seems to self-destruct trying to run on a matrix
  for(i = 0; i < N; i++){
    for(c = 0; c < C; c++){
      sd_ret(c, i) = sqrt((post_sums_sq(c, i) - (num_samp * post_sums(c, i) * post_sums(c, i))) / (num_samp - 1.0));
    }
  }
  sd.push_back(sd_ret);

  mean.names() = CharacterVector::create("pi", "PofZ");
  sd.names() = CharacterVector::create("pi", "PofZ");

  ret = List::create(mean, sd, trace);
  ret.names() = CharacterVector::create("mean", "sd", "trace");
  return(ret);
}
