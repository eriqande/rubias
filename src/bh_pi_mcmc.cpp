#include <Rcpp.h>
using namespace Rcpp;
#include "utilities.h"
#include "rcpp_sampling.h"

//' MCMC from a hierarchical GSI model for rho, pi, and the individual posterior probabilities,
//' with misassignment-scaling
//'
//' Using a matrix of scaled likelihoods, this function samples the posteriors and values of rho,
//' then samples values of omega scaled by their corresponding rho and the inverse of their
//' average rate of correct assignment.  It returns the output in a list.
//'
//' @keywords internal
//'
//' @param SL  matrix of the scaled likelihoods.  This is should have values for each
//' individual in a column (going down in the rows are values for different populations).
//' @param Omega_init  Starting value for the omega (collection mixing proportion) vector.
//' @param Rho_init Starting value for the rho (reporting unit mixing proportion) vector.
//' @param lambda_rho the prior to be added to the reporting unit allocations, in order to
//' generate pseudo-count Dirichlet parameters for the simulation of a new rho vector
//' @param lambda_omega the prior to be added to the collection allocations, in order to
//' generate pseudo-count Dirichlet parameters for the simulation of a new omega vector
//' @param reps total number of reps (sweeps) to do.
//' @param burn_in how many reps to discard in the beginning when doing the mean calculation.
//' They will still be returned in the traces if desired
//' @param sample_int_omega the number of reps between samples being taken for omega
//' traces. If 0 no trace samples are taken.
//' @param sample_int_rho the number of reps between samples being taken for rho.
//' If 0 no trace samples are taken.
//' @param sample_int_PofZ the number of reps between samples being taken for the posterior
//' traces of each individual's collection of origin. If 0 no trace samples are taken.
//' @param sample_int_PofR the number of reps between samples being taken for the posterior
//' traces of each individual's reporting unit of origin. If 0 no trace samples are taken.
//' @param RU_starts a vector of length(rho.size()) + 1, where each element delineates
//' the starting index of a reporting unit in RU_vec (last element is total # collections)
//' @param RU_vec a vector of collection indices, grouped by reporting unit, with groups
//' delineated in RU_starts
//' @param coll2correctRU a vector of average rates at which fish from each collection
//' are assigned to itself or to another collection in the same reporting unit;
//' collections should be in the same order as RU_vec.
//'
//' @examples
//' \dontrun{
//' params <- tcf2param_list(alewife, 17)
//' logl <- geno_logL(params)
//' SL <- apply(exp(logl), 2, function(x) x/sum(x))
//' avg_correct <- avg_coll2correctRU(SL, params$coll,params$RU_starts, params$RU_vec)
//' lambda_omega <- rep(1/params$C, params$C)
//' lambda_rho <- rep(1/(length(params$RU_starts)-1), length(params$RU_starts)-1 )
//' test_bh_mcmc <- gsi_mcmc_2(SL,
//'                            lambda_rho,
//'                            lambda_omega,
//'                            lambda_rho,
//'                            lambda_omega,
//'                            10000, 2500, 50, 50, 50, 50,
//'                            params$RU_starts,
//'                            params$RU_vec,
//'                            avg_correct)
//' }
//'
//' @return \code{gsi_mcmc_2} returns a nested list of MCMC results.
//'
//' \code{$mean} records the mean
//' sampled values for rho and omega in vectors, as well as a matrix of the posterior probability of
//' assignment for every in individual (column) to a collection (PofZ, rows) or reporting unit (PofR, rows)
//'
//' \code{$sd} records the standard deviations for the same values. Sampling for both \code{sd}
//' and \code{mean} are only begun after the burn-in period.
//'
//' \code{$trace} is a list, with each element being a list of samples for the relevant variable
//' (rho, omega, PofZ, PofR) taken at the chosen sampling interval. If the sampling interval for
//' any parameter = 0, that list is empty.
//'
//' @export
// [[Rcpp::export]]
List gsi_mcmc_2(NumericMatrix SL, NumericVector Rho_init, NumericVector Omega_init, NumericVector lambda_rho, NumericVector lambda_omega, int reps, int burn_in, int sample_int_omega, int sample_int_rho, int sample_int_PofZ, int sample_int_PofR, IntegerVector RU_starts, IntegerVector RU_vec, NumericVector coll2correctRU) {
  List omega_list;
  List rho_list;
  List PofZ_list;
  List PofR_list;
  List trace, mean, sd, ret;
  IntegerVector allocs(SL.ncol());
  IntegerVector counts(SL.nrow());
  IntegerVector ru_counts(Rho_init.size());
  NumericVector omega = clone(Omega_init);
  NumericVector omega_sums(Omega_init.size());
  NumericVector omega_sums_sq(Omega_init.size());
  NumericVector rho = clone(Rho_init);
  NumericVector rho_sums(Rho_init.size());
  NumericVector rho_sums_sq(Rho_init.size());
  NumericMatrix posts = clone(SL);
  NumericMatrix post_sums(SL.nrow(), SL.ncol());
  NumericMatrix post_sums_sq(SL.nrow(), SL.ncol());
  NumericMatrix Z_sd_ret(SL.nrow(), SL.ncol());
  NumericMatrix posts_R(rho.size(), SL.ncol());
  NumericMatrix post_R_sums(rho.size(), SL.ncol());
  NumericMatrix post_R_sums_sq(rho.size(), SL.ncol());
  NumericMatrix R_sd_ret(rho.size(), SL.ncol());

  int R = SL.nrow();
  int C = SL.ncol();
  int K = posts_R.nrow();
  int i, r, c, f, coll, ru;
  double sum, tmp;
  int num_samp = reps - burn_in;

  if(num_samp <= 1) stop("reps - burn_in <= 1");


  for(i = 0; i < reps; i++) {
    // store rho values
    if( (sample_int_rho > 0) && (i % sample_int_rho == 0) ) {
      rho_list.push_back(rho);
    }
    if(i >= burn_in) {
      rho_sums += rho;
      rho_sums_sq += rho * rho;
    }
    // store omega value
    if( (sample_int_omega > 0) && (i % sample_int_omega == 0) ) {
      omega_list.push_back(omega);
    }
    if(i >= burn_in) {
      omega_sums += omega;
      omega_sums_sq += omega * omega;
    }
    // normalize the scaled likelihoods into posteriors
    for(c = 0; c < C; c++) {
      sum = 0.0;
      for(r = 0; r < R; r++) {
        tmp = SL(r, c) * omega[r];
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

    //calculate PofR values, if past burn_in
    if( i >= burn_in ) {
      // initialize posts_R to 0
      for(c = 0; c < C; c++) {
        for(ru = 0; ru < K; ru++) {
          posts_R(ru,c) = 0;
        }
      }
      for(c = 0; c < C; c++) {
        sum = 0.0;
        for(ru = 0; ru < K; ru++) {
          for(coll = RU_starts[ru]; coll < RU_starts[ru + 1]; coll++) {
            posts_R(ru,c) += posts(RU_vec[coll]-1, c);
          }
        }
        for(ru = 0; ru < K; ru++) {
          post_R_sums(ru, c) += posts_R(ru,c);
          post_R_sums_sq(ru, c) += posts_R(ru, c) * posts_R(ru, c);
        }

      }
      if( (sample_int_PofR > 0) && (i % sample_int_PofR == 0) ) {
        PofR_list.push_back(posts_R);
      }
    }

    // allocate individuals to populations
    allocs = samp_from_mat(posts);
    // tabulate the counts of individuals per population
    for(r = 0; r < R; r++) {
      counts[r] = 0;
    }
    for(f = 0; f < C; f++) {
      counts[allocs[f] - 1] += 1;
    }
    // condense to reporting units
    for(ru = 0; ru < K; ru++) {
      ru_counts[ru] = 0;
    }
    for(ru = 0; ru < K; ru++){
      for(coll = RU_starts[ru]; coll < RU_starts[ru + 1]; coll++) {
        ru_counts[ru] += counts[RU_vec[coll]-1];
      }
    }
    // simulate new rho values
    rho = dirch_from_counts(ru_counts, lambda_rho);
    // simulate omega parameters for all collections in a reporting unit, and weight by rho
    sum = 0.0;
    for(ru = 0; ru < K; ru++){
      int C = RU_starts[ru+1] - RU_starts[ru];
      IntegerVector tmp_counts(C);
      NumericVector tmp_lambdas(C);
      NumericVector tmp_oms(C);
      coll = 0;
      for(c = RU_starts[ru]; c < RU_starts[ru + 1]; c++) {
        tmp_counts[coll] = counts[RU_vec[c] - 1];
        tmp_lambdas[coll] = lambda_omega[RU_vec[c] - 1];
        coll += 1;
      }
      tmp_oms = dirch_from_counts(tmp_counts, tmp_lambdas);
      coll = 0;
      for(c = RU_starts[ru]; c < RU_starts[ru + 1]; c++) {
        omega[RU_vec[c] - 1] = tmp_oms[coll] * rho[ru] * 1/coll2correctRU[c];
        sum += omega[RU_vec[c] - 1];
        coll += 1;
      }
    }
    for(coll = 0; coll < R; coll ++) {
      omega[coll] /= sum;
    }
  }

  // put the traces in there if there are any
  trace = List::create(rho_list, PofR_list, omega_list, PofZ_list);
  trace.names() = CharacterVector::create("rho","PofR","omega", "PofZ");

  rho_sums = rho_sums / num_samp;
  post_R_sums = post_R_sums / num_samp;
  post_sums = post_sums / num_samp;
  omega_sums  = omega_sums / num_samp;

  mean.push_back(rho_sums);
  sd.push_back(sqrt((rho_sums_sq - (num_samp * rho_sums * rho_sums)) / (num_samp - 1.0)));

  mean.push_back(post_R_sums);
  for(c = 0; c < C; c++){
    for(ru = 0; ru < K; ru++){
      R_sd_ret(ru,c) = sqrt((post_R_sums_sq(ru,c) - (num_samp * post_R_sums(ru,c) * post_R_sums(ru,c))) / (num_samp - 1.0));
    }
  }
  sd.push_back(R_sd_ret);

  mean.push_back(omega_sums);
  sd.push_back(sqrt((omega_sums_sq - (num_samp * omega_sums * omega_sums)) / (num_samp - 1.0)));

  mean.push_back(post_sums);
  for(c = 0; c < C; c++){
    for(r = 0; r < R; r++){
      Z_sd_ret(r,c) = sqrt((post_sums_sq(r,c) - (num_samp * post_sums(r,c) * post_sums(r,c))) / (num_samp - 1.0));
    }
  }
  sd.push_back(Z_sd_ret);

  mean.names() = CharacterVector::create("rho", "PofR", "omega", "PofZ");
  sd.names() = CharacterVector::create("rho", "PofR", "omega", "PofZ");

  ret = List::create(mean, sd, trace);
  ret.names() = CharacterVector::create("mean", "sd", "trace");
  return(ret);
}







////////////////////////////
//' MCMC from a hierarchical GSI model for rho, pi, and the individual posterior probabilities,
//' with misassignment-scaling
//'
//' This is a mild re-write.  Eric wanted to remove the coll2correctRU scaling.
//'
//' Using a matrix of scaled likelihoods, this function samples the posteriors and values of rho,
//' then samples values of omega scaled by their corresponding rho and the inverse of their
//' average rate of correct assignment.  It returns the output in a list.
//'
//'
//' @param SL  matrix of the scaled likelihoods.  This is should have values for each
//' individual in a column (going down in the rows are values for different populations).
//' @param Omega_init  Starting value for the omega (collection mixing proportion) vector.
//' @param Rho_init Starting value for the rho (reporting unit mixing proportion) vector.
//' @param lambda_rho the prior to be added to the reporting unit allocations, in order to
//' generate pseudo-count Dirichlet parameters for the simulation of a new rho vector
//' @param lambda_omega the prior to be added to the collection allocations, in order to
//' generate pseudo-count Dirichlet parameters for the simulation of a new omega vector
//' @param reps total number of reps (sweeps) to do.
//' @param burn_in how many reps to discard in the beginning when doing the mean calculation.
//' They will still be returned in the traces if desired
//' @param sample_int_omega the number of reps between samples being taken for omega
//' traces. If 0 no trace samples are taken.
//' @param sample_int_rho the number of reps between samples being taken for rho.
//' If 0 no trace samples are taken.
//' @param sample_int_PofZ the number of reps between samples being taken for the posterior
//' traces of each individual's collection of origin. If 0 no trace samples are taken.
//' @param sample_int_PofR the number of reps between samples being taken for the posterior
//' traces of each individual's reporting unit of origin. If 0 no trace samples are taken.
//' @param RU_starts a vector of length(rho.size()) + 1, where each element delineates
//' the starting index of a reporting unit in RU_vec (last element is total # collections)
//' @param RU_vec a vector of collection indices, grouped by reporting unit, with groups
//' delineated in RU_starts
//'
//' @examples
//' \dontrun{
//' params <- tcf2param_list(alewife, 15)
//' logl <- geno_logL(params)
//' SL <- apply(exp(logl), 2, function(x) x/sum(x))
//' lambda_omega <- rep(1/params$C, params$C)
//' lambda_rho <- rep(1/(length(params$RU_starts)-1), length(params$RU_starts)-1 )
//' test_bh_mcmc <- gsi_mcmc_bh(SL,
//'                             lambda_rho,
//'                             lambda_omega,
//'                             lambda_rho,
//'                             lambda_omega,
//'                             10000, 2500, 50, 50, 50, 50,
//'                             params$RU_starts,
//'                             params$RU_vec)
//' }
//'
//' @keywords internal
//' @return \code{gsi_mcmc_2} returns a nested list of MCMC results.
//'
//' \code{$mean} records the mean
//' sampled values for rho and omega in vectors, as well as a matrix of the posterior probability of
//' assignment for every in individual (column) to a collection (PofZ, rows) or reporting unit (PofR, rows)
//'
//' \code{$sd} records the standard deviations for the same values. Sampling for both \code{sd}
//' and \code{mean} are only begun after the burn-in period.
//'
//' \code{$trace} is a list, with each element being a list of samples for the relevant variable
//' (rho, omega, PofZ, PofR) taken at the chosen sampling interval. If the sampling interval for
//' any parameter = 0, that list is empty.
//'
//' @export
// [[Rcpp::export]]
List gsi_mcmc_bh(NumericMatrix SL, NumericVector Rho_init, NumericVector Omega_init, NumericVector lambda_rho, NumericVector lambda_omega, int reps, int burn_in, int sample_int_omega, int sample_int_rho, int sample_int_PofZ, int sample_int_PofR, IntegerVector RU_starts, IntegerVector RU_vec) {
  List omega_list;
  List rho_list;
  List PofZ_list;
  List PofR_list;
  List trace, mean, sd, ret;
  IntegerVector allocs(SL.ncol());
  IntegerVector counts(SL.nrow());
  IntegerVector ru_counts(Rho_init.size());
  NumericVector omega = clone(Omega_init);
  NumericVector omega_sums(Omega_init.size());
  NumericVector omega_sums_sq(Omega_init.size());
  NumericVector rho = clone(Rho_init);
  NumericVector rho_sums(Rho_init.size());
  NumericVector rho_sums_sq(Rho_init.size());
  NumericMatrix posts = clone(SL);
  NumericMatrix post_sums(SL.nrow(), SL.ncol());
  NumericMatrix post_sums_sq(SL.nrow(), SL.ncol());
  NumericMatrix Z_sd_ret(SL.nrow(), SL.ncol());
  NumericMatrix posts_R(rho.size(), SL.ncol());
  NumericMatrix post_R_sums(rho.size(), SL.ncol());
  NumericMatrix post_R_sums_sq(rho.size(), SL.ncol());
  NumericMatrix R_sd_ret(rho.size(), SL.ncol());

  int R = SL.nrow();
  int C = SL.ncol();
  int K = posts_R.nrow();
  int i, r, c, f, coll, ru;
  double sum, tmp;
  int num_samp = reps - burn_in;

  if(num_samp <= 1) stop("reps - burn_in <= 1");


  for(i = 0; i < reps; i++) {
    // store rho values
    if( (sample_int_rho > 0) && (i % sample_int_rho == 0) ) {
      rho_list.push_back(rho);
    }
    if(i >= burn_in) {
      rho_sums += rho;
      rho_sums_sq += rho * rho;
    }
    // store omega value
    if( (sample_int_omega > 0) && (i % sample_int_omega == 0) ) {
      omega_list.push_back(omega);
    }
    if(i >= burn_in) {
      omega_sums += omega;
      omega_sums_sq += omega * omega;
    }
    // normalize the scaled likelihoods into posteriors
    for(c = 0; c < C; c++) {
      sum = 0.0;
      for(r = 0; r < R; r++) {
        tmp = SL(r, c) * omega[r];
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

    //calculate PofR values, if past burn_in
    if( i >= burn_in ) {
      // initialize posts_R to 0
      for(c = 0; c < C; c++) {
        for(ru = 0; ru < K; ru++) {
          posts_R(ru,c) = 0;
        }
      }
      for(c = 0; c < C; c++) {
        sum = 0.0;
        for(ru = 0; ru < K; ru++) {
          for(coll = RU_starts[ru]; coll < RU_starts[ru + 1]; coll++) {
            posts_R(ru,c) += posts(RU_vec[coll]-1, c);
          }
        }
        for(ru = 0; ru < K; ru++) {
          post_R_sums(ru, c) += posts_R(ru,c);
          post_R_sums_sq(ru, c) += posts_R(ru, c) * posts_R(ru, c);
        }

      }
      if( (sample_int_PofR > 0) && (i % sample_int_PofR == 0) ) {
        PofR_list.push_back(posts_R);
      }
    }

    // allocate individuals to populations
    allocs = samp_from_mat(posts);
    // tabulate the counts of individuals per population
    for(r = 0; r < R; r++) {
      counts[r] = 0;
    }
    for(f = 0; f < C; f++) {
      counts[allocs[f] - 1] += 1;
    }
    // condense to reporting units
    for(ru = 0; ru < K; ru++) {
      ru_counts[ru] = 0;
    }
    for(ru = 0; ru < K; ru++){
      for(coll = RU_starts[ru]; coll < RU_starts[ru + 1]; coll++) {
        ru_counts[ru] += counts[RU_vec[coll]-1];
      }
    }
    // simulate new rho values
    rho = dirch_from_counts(ru_counts, lambda_rho);
    // simulate omega parameters for all collections in a reporting unit, and weight by rho
    sum = 0.0;
    for(ru = 0; ru < K; ru++){
      int C = RU_starts[ru+1] - RU_starts[ru];
      IntegerVector tmp_counts(C);
      NumericVector tmp_lambdas(C);
      NumericVector tmp_oms(C);
      coll = 0;
      for(c = RU_starts[ru]; c < RU_starts[ru + 1]; c++) {
        tmp_counts[coll] = counts[RU_vec[c] - 1];
        tmp_lambdas[coll] = lambda_omega[RU_vec[c] - 1];
        coll += 1;
      }
      tmp_oms = dirch_from_counts(tmp_counts, tmp_lambdas);
      coll = 0;
      for(c = RU_starts[ru]; c < RU_starts[ru + 1]; c++) {
        omega[RU_vec[c] - 1] = tmp_oms[coll] * rho[ru];
        sum += omega[RU_vec[c] - 1];
        coll += 1;
      }
    }
    for(coll = 0; coll < R; coll ++) {
      omega[coll] /= sum;
    }
  }

  // put the traces in there if there are any
  trace = List::create(rho_list, PofR_list, omega_list, PofZ_list);
  trace.names() = CharacterVector::create("rho","PofR","omega", "PofZ");

  rho_sums = rho_sums / num_samp;
  post_R_sums = post_R_sums / num_samp;
  post_sums = post_sums / num_samp;
  omega_sums  = omega_sums / num_samp;

  mean.push_back(rho_sums);
  sd.push_back(sqrt((rho_sums_sq - (num_samp * rho_sums * rho_sums)) / (num_samp - 1.0)));

  mean.push_back(post_R_sums);
  for(c = 0; c < C; c++){
    for(ru = 0; ru < K; ru++){
      R_sd_ret(ru,c) = sqrt((post_R_sums_sq(ru,c) - (num_samp * post_R_sums(ru,c) * post_R_sums(ru,c))) / (num_samp - 1.0));
    }
  }
  sd.push_back(R_sd_ret);

  mean.push_back(omega_sums);
  sd.push_back(sqrt((omega_sums_sq - (num_samp * omega_sums * omega_sums)) / (num_samp - 1.0)));

  mean.push_back(post_sums);
  for(c = 0; c < C; c++){
    for(r = 0; r < R; r++){
      Z_sd_ret(r,c) = sqrt((post_sums_sq(r,c) - (num_samp * post_sums(r,c) * post_sums(r,c))) / (num_samp - 1.0));
    }
  }
  sd.push_back(Z_sd_ret);

  mean.names() = CharacterVector::create("rho", "PofR", "omega", "PofZ");
  sd.names() = CharacterVector::create("rho", "PofR", "omega", "PofZ");

  ret = List::create(mean, sd, trace);
  ret.names() = CharacterVector::create("mean", "sd", "trace");
  return(ret);
}




