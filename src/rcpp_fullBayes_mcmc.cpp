#include <Rcpp.h>
#include <RcppParallel.h>
#include "utilities.h"
#include "macros.h"
#include "rcpp_sampling.h"
using namespace Rcpp;
using namespace RcppParallel;


//' MCMC from the fully Bayesian GSI model for pi and the individual posterior probabilities
//'
//' Given a list of key parameters from a genetic dataset, this function samples values of pi
//' and the posteriors for all the individuals. Each MCMC iteration includes a recalculation
//' of the scaled genotype likelihood matrix, with baseline allele frequencies updated
//' based on the previous iteration's allocations. It returns the output in a list.
//' @keywords internal
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//' @param Pi_init  Starting value for the pi (collection mixture proportion) vector.
//' @param lambda the prior to be added to the collection allocations, in order to
//' generate pseudo-count Dirichlet parameters for the simulation of a new pi vector
//' @param reps total number of reps (sweeps) to do.
//' @param burn_in how many reps to discard in the beginning when doing the mean calculation. They will still be
//' returned in the traces if desired
//' @param sample_int_Pi the number of reps between samples being taken for Pi traces.  If 0 no trace samples are taken
//' @param sample_int_PofZ the number of reps between samples being taken for the traces of posterior of each individual's origin. If 0
//' no trace samples are taken.
//'  @param sample_total_catch integer. Set it to 1 if you want to sample the total-stock specific catch. If so,
//' then you have to set `total_catch_vals` appropriately.
//' @param total_catch_vals an integer vector of length reps that holds the total catch.  It is a vector to allow
//' for this to be a sample from the posterior for the total catch.
//' @param variable_prob_is_catch integer. Set to 1 if samples have different probabilities of being
//' considered catch.  If it is 1, then \code{prob_is_catch_vec} must be provided.
//' @param prob_is_catch_vec NumericVector of probabilities that each individual in the sample should
//' be considered catch.
//'
//' @return \code{gsi_mcmc_fb} returns a list of three. \code{$mean} lists the posterior
//' means for collection proportions \code{pi}, for the individual posterior
//' probabilities of assignment \code{PofZ}, and for the allele frequencies \code{theta}.
//' \code{$sd} returns the posterior standard deviations for the same values.
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
//' lambda <- rep(1/params$C, params$C)
//' # use very short run and burn in so it doesn't take too long
//' # when checking on CRAN
//' mcmc <- gsi_mcmc_fb(params, lambda, lambda, 20, 5, 4, 4)
//' @export
// [[Rcpp::export]]
List gsi_mcmc_fb(
    List par_list,
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

  // Code from geno_logl, creating C X N loglikelihood matrix
  int r, i, c, l, a1, a2;
  double sum, tmp;
  int N = as<int>(par_list["N"]);
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  IntegerVector coll = as<IntegerVector>(par_list["coll"]);
  IntegerVector PLOID = as<IntegerVector>(par_list["ploidies"]);
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  NumericVector theta(DP.size());
  NumericMatrix logl(C, N);

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
  NumericVector theta_sums(DP.size());
  NumericVector theta_sums_sq(DP.size());
  NumericMatrix posts = clone(logl);
  NumericMatrix post_sums(logl.nrow(), logl.ncol());
  NumericMatrix post_sums_sq(logl.nrow(), logl.ncol());
  NumericMatrix sd_ret(logl.nrow(), logl.ncol());
  IntegerVector allocs(logl.ncol());
  NumericVector DP_temp = clone(DP);
  NumericVector sum_DP_temp = clone(sum_DP);

  int num_samp = reps - burn_in;
  if(num_samp <= 1) stop("reps - burn_in <= 1");


  // Creating Worker object for RcppParallel
  struct GenoLike : public Worker
  {
    // source datasets
    const int N;
    const int C;
    const int L;
    const RVector<int> A;
    const RVector<int> CA;
    const RVector<int> coll;
    const RVector<int> PLOID;
    const RVector<int> I;
    const RVector<double> theta;

    // destination matrix
    RMatrix<double> logl;

    GenoLike(const int N, const int C, const int L, const IntegerVector A, const IntegerVector CA,
             const IntegerVector coll, const IntegerVector PLOID,
             const IntegerVector I, const NumericVector theta, NumericMatrix logl)
      : N(N), C(C), L(L), A(A), CA(CA), coll(coll), PLOID(PLOID), I(I), theta(theta), logl(logl) {}

    void operator()(std::size_t begin, std::size_t end) {
      int c, l;
      double sum, y1, y2, gp, colmean, colsum;
      for (std::size_t i = begin; i < end; i++) {
        colsum = 0.0;
        for(c = 0; c < C; c++) { // cycle over collections
          sum = 0.0;
          for(l = 0; l < L; l++) {  // cycle over loci
            int a1 = I[I_dx(l, i, 0, 2, N)] - 1;
            int a2 = I[I_dx(l, i, 1, 2, N)] - 1;
            if(PLOID[l] == 1) {
              if(a1 < 0) {gp = 1.0;} else {gp = theta[D_dx(l, c, a1, L, C, A, CA)];}
            }
            else {
              if(a1 < 0 || a2 < 0) {gp = 1.0;} else {
                y1 = theta[D_dx(l, c, a1, L, C, A, CA)];
                y2 = theta[D_dx(l, c, a2, L, C, A, CA)];
                gp = y1 * y2 * (1 + (a1 == a2));
              }}
            sum += log(gp);
          }
          logl(c, i) = sum;    // Trim down to just a c-long vector?
          colsum += sum; // sum across collections for column mean calculation
        }

        // take column means, then sweep out from each logl to prevent underflow
        colmean = colsum/C;
        for(c = 0; c < C; c++) {
          logl(c, i) = logl(c, i) - colmean;
        }
      }
    }
  };

  //begin MCMC cycling
  for(r = 0; r < reps; r++) {

    checkUserInterrupt();
    // store pi value
    if( (sample_int_Pi > 0) && (r % sample_int_Pi == 0) ) {
      pi_list.push_back(pi);
    }
    if(r >= burn_in) {
      pi_sums += pi;
      pi_sums_sq += pi * pi;
    }

    // simulate allele frequency values based on Dirichlet parameter vector
    for(c = 0; c < C; c++) { // cycle over collections
      for(l = 0; l < L; l++) { // cycle over loci
      sum = 0.0;
        for(a1 = 0; a1 < A[l]; a1++) { // cycle over alleles within the locus
          tmp = rgammadouble(1L, DP_temp[D_dx(l, c, a1, L, C, A, CA)], 1.0);
          theta[D_dx(l, c, a1, L, C, A, CA)] = tmp;
          sum += tmp;
        }
        for(a1 = 0; a1 < A[l]; a1++) { // cycle again to normalize to 1
          theta[D_dx(l, c, a1, L, C, A, CA)] /= sum;
        }
      }
    }

    // store theta value
    if(r >= burn_in) {
      theta_sums += theta;
      theta_sums_sq += theta * theta;
    }

    // genotype likelihood calculations

    GenoLike genoLike(N, C, L, A, CA, coll, PLOID, I, theta, logl);
    parallelFor(0, logl.ncol(), genoLike);


    // convert to likelihood and calculate the posterior in one go
    for(i = 0; i < N; i++) {
      sum = 0.0;
      for(c = 0; c < C; c++) {
        tmp = exp(logl(c,i)) * pi[c];
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


    // Copy allocs into Zeds for stock-specific total catch stuff
    IntegerVector Zeds = clone(allocs);


    pi = dirch_from_allocations(allocs, lambda);

    // compute a new Dirichlet Parameter Vector based on the allocations

    std::copy(DP.begin(), DP.end(), DP_temp.begin());
    std::copy(sum_DP.begin(), sum_DP.end(), sum_DP_temp.begin());

    for(i = 0; i < N; i++) {
      c = allocs[i] - 1;
      for(l = 0; l < L; l++) {
        a1 = I[I_dx(l, i, 0, 2, N)] - 1;
        a2 = I[I_dx(l, i, 1, 2, N)] - 1;


        if(PLOID[l] == 1) {
          if(a1 >= 0) {
            DP_temp[D_dx(l, c, a1, L, C, A, CA)] += 1;
            sum_DP_temp[SD_dx(l, c, C)] += 1;
          }


        } else {
          if(a1 >= 0 && a2 >= 0) {
            DP_temp[D_dx(l, c, a1, L, C, A, CA)] += 1;
            DP_temp[D_dx(l, c, a2, L, C, A, CA)] += 1;
            sum_DP_temp[SD_dx(l, c, C)] += 2;
          }


        }
      }
    } // Done with computing a new Dirichlet Parameter Vector


    // Now, deal with all the stock-specific total catch stuff
    if(sample_total_catch != 0) {

      int NumThatAreCatch = N;  // by default, everyone is catch, N = number of samples

      // now, if we have variable_prob_is_catch we need to simulate whether
      // these individuals are part of the catch or not.  If they are not
      // part of the catch, we assign them -1's. This also updates NumThatAreCatch
      if(variable_prob_is_catch == 1) {
        Zeds = turn_non_catch_to_minus_one(Zeds, prob_is_catch_vec, NumThatAreCatch);
      }


      //if(i % 500 == 0)  Rcpp::Rcout << "i= " << i << "  NewZeds= " << Zeds << std::endl;


      // count up the current allocations (the catch contribution of the actual samples)
      IntegerVector CurrentAllocations = tabulate_allocations(Zeds, C);



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

    }  // end of stock-specific catch block



  }  // close loop over r (reps)

  // put the traces in there if there are any
  if(sample_total_catch == 0) {
    trace = List::create(pi_list, PofZ_list);
    trace.names() = CharacterVector::create("pi", "PofZ");
  } else {
    trace = List::create(pi_list, PofZ_list, SSC_list, PPRC_list, CA_list);
    trace.names() = CharacterVector::create("pi", "PofZ", "SSTC", "PPRC", "CA");
  }

  // put the means and standard devs and traces in the return variable
  post_sums = post_sums / num_samp;
  pi_sums  = pi_sums / num_samp;
  theta_sums = theta_sums / num_samp;

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

  mean.push_back(theta_sums);
  sd.push_back(sqrt((theta_sums_sq - (num_samp * theta_sums * theta_sums)) / (num_samp - 1.0)));

  mean.names() = CharacterVector::create("pi", "PofZ", "theta");
  sd.names() = CharacterVector::create("pi", "PofZ", "theta");

  ret = List::create(mean, sd, trace);
  return(ret);
}
