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
//' @param sample_total_catch integer. Set it to 1 if you want to sample the total-stock specific catch. If so,
//' then you have to set `total_catch_vals` appropriately.
//' @param total_catch_vals an integer vector of length reps that holds the total catch.  It is a vector to allow
//' for this to be a sample from the posterior for the total catch.
//' @param variable_prob_is_catch integer. Set to 1 if samples have different probabilities of being
//' considered catch.  If it is 1, then \code{prob_is_catch_vec} must be provided.
//' @param prob_is_catch_vec NumericVector of probabilities that each individual in the sample should
//' be considered catch.
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
List gsi_mcmc_1(
    NumericMatrix SL,
    NumericVector Pi_init,
    NumericVector lambda,
    int reps,
    int burn_in,
    int sample_int_Pi,
    int sample_int_PofZ,
    int sample_total_catch = 0,
    IntegerVector total_catch_vals = IntegerVector::create(-1),  // this makes the default a vector of length 1, holding -1.  Hopefully this will work when the user does not give a value for this...
    int variable_prob_is_catch = 0,
    NumericVector prob_is_catch_vec = NumericVector::create(-1)
) {
  List pi_list;
  List PofZ_list;
  List SSC_list;   // For gathering stock-specific total catch over the iterations
  List PPRC_list;  // (PPRC = posterior-predictive remaining catch) for gathering stock-specific remaining catch over iterations (i.e., the part simulated from the posterior predictive of the non-sampled fish)
  List CA_list;   // (CA = Current Allocations) for gathering a summary of the current allocations over the iterations.
                  // Note: SSC = CA + PPRC.  I just want to be able to return it broken out by CA and PPRC

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

    // allocate individuals to populations and simulate a new pi.
    // This was updated in August 2025 to store the allocations (the Z's are Zeds)
    // of the fish in the mixture.  Doing so gives us access to them, which
    // allows us to sample from the posterior predictive of the total catch.
    IntegerVector Zeds = samp_from_mat(posts);

    //if(i % 500 == 0)  Rcpp::Rcout << "i= " << i << "     Zeds= " << Zeds << std::endl;

    // finally, here we update pi
    pi = dirch_from_allocations(Zeds, lambda);


    // If we are sampling from the posterior predictive of remaining catch
    // we do it here
    if(sample_total_catch != 0) {

      int NumThatAreCatch = C;  // by default, everyone is catch, C = number of columns in SL matrix = number of samples.

      // now, if we have variable_prob_is_catch we need to simulate whether
      // these individuals are part of the catch or not.  If they are not
      // part of the catch, we assign them -1's. This also updates NumThatAreCatch
      if(variable_prob_is_catch == 1) {
        Zeds = turn_non_catch_to_minus_one(Zeds, prob_is_catch_vec, NumThatAreCatch);
      }


      //if(i % 500 == 0)  Rcpp::Rcout << "i= " << i << "  NewZeds= " << Zeds << std::endl;


      // count up the current allocations (the catch contribution of the actual samples)
      IntegerVector CurrentAllocations = tabulate_allocations(Zeds, R);



      int remaining_catch = total_catch_vals(i) - NumThatAreCatch;
      if(remaining_catch < 0) remaining_catch = 0;

      unsigned int ncell = pi.length();
      unsigned int rem_catch_uns = remaining_catch;

      //if(i % 500 == 0)  Rcpp::Rcout << "i= " << i << "  NumThatAreCatch= " << NumThatAreCatch << "   rem_catch_uns= " << rem_catch_uns << std::endl;

      // simulate the posterior predictive of the "remaining catch", PPRC
      IntegerVector PPRC = rmultinom_1(rem_catch_uns, pi, ncell);

      //if(i % 500 == 0) Rcpp::Rcout << "i= " << i << "  PPRC= " << PPRC << std::endl;


      //if(i % 500 == 0) Rcpp::Rcout << "i= " << i << "  CUAL= " << CurrentAllocations << std::endl;

      SSC_list.push_back(PPRC + CurrentAllocations);

      // Hey Eric!  Return the SSC_list as done here, but also you are going to
      // want to return the Remaining and the CurrentAllocations separately...
      PPRC_list.push_back(PPRC);
      CA_list.push_back(CurrentAllocations);

    }


  }

  // put the traces in there if there are any
  if(sample_total_catch == 0) {
    trace = List::create(pi_list, PofZ_list);
    trace.names() = CharacterVector::create("pi", "PofZ");
  } else {
    trace = List::create(pi_list, PofZ_list, SSC_list, PPRC_list, CA_list);
    trace.names() = CharacterVector::create("pi", "PofZ", "SSTC", "PPRC", "CA");
  }
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
