#include <Rcpp.h>
#include "macros.h"
#include "utilities.h"
using namespace Rcpp;

//' Simulate genotype log-likelihoods from a population by gene copy
//'
//' Takes a list of parameters from a genetic dataset, and returns a genotype log-likelihood
//' matrix for individuals simulated by gene copy from the specified collections
//'
//' In simulation by gene copy, the genotype at a locus for any individual is the result
//' of two random draws from the allele count matrix of that locus. Draws within an individual
//' are performed without replacement, but allele counts are replaced between individuals.
//'
//' @keywords internal
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//' @param sim_colls a vector of indices for the collections desired for simulation;
//' each element of the list corresponds to an individual
//'
//'
//' @return \code{gprob_sim} returns a matrix of the summed log-likelihoods
//' for all loci of a simulated population mixture; columns represent individuals,
//' with each row containing their log-likelihood of belonging to the collection
//' of the same index, given the selection of two independent gene copies from the
//' desired collection of origin's reference allele frequencies
//'
//' @examples
//' example(tcf2param_list)
//' sim_colls <- sample(ale_par_list$C, 1070, replace = T)
//' ale_sim_gprobs_gc <- gprob_sim_gc(ale_par_list, sim_colls)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix gprob_sim_gc(List par_list, IntegerVector sim_colls) {
  int i, c, l , a, S, S2;
  int N = sim_colls.size();
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);
  int LOO;
  IntegerVector a1_vec(L);
  IntegerVector a2_vec(L);
  IntegerVector AC = as<IntegerVector>(par_list["AC"]);
  IntegerVector sum_AC = as<IntegerVector>(par_list["sum_AC"]);
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  double cumul, rando, sum, gp;
  NumericMatrix out(C,N);


  for(i = 0; i < N; i++) {
    for(l = 0; l < L; l++) {
      S = sum_AC[SD_dx(l, sim_colls[i]-1, C)];
      S2 = S - 1;
      if(S == 0 || S == 1) {
        a1_vec[l] = -1;
        a2_vec[l] = -1;
      } else {
        rando = runif(1)[0] * S;
        cumul = 0.0;
        for(a = 0; a < A[l]; a++) {
          cumul += AC[D_dx(l, sim_colls[i]-1, a, L, C, A, CA)];
          a1_vec[l] = a;
          if (cumul >= rando) {
            break;
          }
        }

        rando = runif(1)[0] * S2;
        cumul = 0.0;
        for(a = 0; a < A[l]; a++) {
          cumul += AC[D_dx(l, sim_colls[i]-1, a, L, C, A, CA)] - (a1_vec[l] == a);
          a2_vec[l] = a;
          if (cumul >= rando) {
            break;
          }
        }
      }
    }
      for(c = 0; c < C; c++) {
        sum = 0.0;
        LOO = c == (sim_colls[i]-1);
        for(l = 0; l < L; l++) {
          GPROB_DIP_FROM_SIM(a1_vec[l], a2_vec[l], l, c, gp);
          sum += log(gp);
        }
        out(c, i) = sum;
      }
  }
  return(out);
}

//' Simulate genotype log-likelihoods from a population by individual
//'
//' Takes a list of parameters from a genetic dataset, and returns a genotype log-likelihood
//' matrix for individuals simulated by individual from the specified collections
//'
//' In simulation by individual, the genotype for any simulated individual is the
//' result of a single random draw from the genotypes of all individuals in the collection.
//' Gene copies and loci are therefore not independent.
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//' @param sim_colls a vector of indices for the collections desired for simulation;
//' each element of the list corresponds to an individual
//'
//'
//' @return \code{gprob_sim} returns a matrix of the summed log-likelihoods
//' for all loci of a simulated population mixture; columns represent individuals,
//' with each row containing their log-likelihood of belonging to the collection
//' of the same index, given the selection of an individual's genotype from the
//' reference collection of interest. Selection at the locus and gene copy level
//' are not independent, and missing data is included in selection.
//'
//' @keywords internal
//'
//' @examples
//' example(tcf2param_list)
//' sim_colls <- sample(ale_par_list$C, 1070, replace = T)
//' ale_sim_gprobs_ind <- gprob_sim_ind(ale_par_list, sim_colls)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix gprob_sim_ind(List par_list, IntegerVector sim_colls) {
  int i, c, l, k, n, ind, count;
  int N = sim_colls.size();
  int K = as<int>(par_list["N"]);
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);
  int LOO;
  IntegerVector a1_vec(L);
  IntegerVector a2_vec(L);
  IntegerVector coll_N = as<IntegerVector>(par_list["coll_N"]);
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  IntegerVector I_coll = as<IntegerVector>(par_list["coll"]);
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  double rando = 0.0, sum, gp;
  NumericMatrix out(C,N);


  for(i = 0; i < N; i++) {
    n = coll_N[sim_colls[i] - 1];
    ind = randint(n);
    count = 0;
    for(k = 0; k < K; k++) {           // Not a terribly efficient way to get the (ind)th individual of a collection, but best way without knowing I is organized by collection, and without a new indexing macro
      if(I_coll[k] == sim_colls[i]) {
        count += 1;
      }
      if(ind == count) {
        break;
      }
    }
    for(l = 0; l < L; l++) {
        a1_vec[l] = I[I_dx(l, k, 0, 2, K)] - 1;
        a2_vec[l] = I[I_dx(l, k, 1, 2, K)] - 1;
      }
    for(c = 0; c < C; c++) {
      sum = 0.0;
      LOO = c == (sim_colls[i]-1);
      for(l = 0; l < L; l++) {
        GPROB_DIP_FROM_SIM(a1_vec[l], a2_vec[l], l, c, gp);
        sum += log(gp);
      }
      out(c, i) = sum;
    }
  }
  return(out);
}

//' Simulate genotypes by gene copy, with missing data from chosen individuals
//'
//' Takes a list of parameters from a genetic dataset, and returns a genotype log-likelihood
//' matrix for individuals simulated by gene copy from the specified collections, with
//' genotypes masked by missing data patterns from reference individuals
//'
//' In simulation by gene copy, the genotype at a locus for any individual is the result
//' of two random draws from the allele count matrix of that locus. Draws within an individual
//' are performed without replacement, but allele counts are replaced between individuals.
//' If the data at a particular locus is missing for individual i in \code{sim_missing},
//' this data will also be missing in simulated individual i for the
//' log-likelihood calculation.
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//' @param sim_colls a vector; element i specifies the collection from which to sample
//' the genotypes for individual i
//' @param sim_missing a vector; element i specifies the index for the individual in
//' params$I whose missing data should be copied for individual i
//'
//' @keywords internal
//' @examples
//' # If one wanted to simulate the missing data patterns
//' of a troublesome mixture dataset, one would run tcf2param_list,
//' selecting samp_type = "mixture", then draw sim_miss from
//' the mixture individual genotype list
//'
//' drawn <- mixture_draw(alewife, rhos = c(1/3, 1/3, 1/3),N = 100)
//' ref <- drawn$reference
//' mix <- drawn$mix
//' params <- tcf2param_list(rbind(ref,mix), 15, samp_type = "mixture")
//' sim_colls <- sample(params$C, 1070, replace = T)
//' sim_miss <- sample(length(params$coll), 1070, replace = T)
//' ale_sim_gprobs_miss <- gprob_sim_gc_missing(params, sim_colls, sim_miss)
//' @export
// [[Rcpp::export]]
NumericMatrix gprob_sim_gc_missing(List par_list, IntegerVector sim_colls, IntegerVector sim_missing) {
  int i, c, l , a, S, S2;
  int N = as<int>(par_list["N"]);
  int K = sim_colls.size();
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);

  int LOO;
  IntegerVector a1_vec(L);
  IntegerVector a2_vec(L);
  IntegerVector AC = as<IntegerVector>(par_list["AC"]);
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  IntegerVector sum_AC = as<IntegerVector>(par_list["sum_AC"]);
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  double cumul, rando, sum, gp;
  NumericMatrix out(C,K);


  for(i = 0; i < K; i++) {
    for(l = 0; l < L; l++) {
      S = sum_AC[SD_dx(l, sim_colls[i]-1, C)];
      S2 = S - 1;
      if(S == 0 || S == 1) {
        a1_vec[l] = -1;
        a2_vec[l] = -1;
      } else {
        if(I[I_dx(l, sim_missing[i]-1, 0, 2, N)] == 0) {
          a1_vec[l] = -1;
        } else {
          rando = runif(1)[0] * S;
          cumul = 0.0;
          for(a = 0; a < A[l]; a++) {
            cumul += AC[D_dx(l, sim_colls[i]-1, a, L, C, A, CA)];
            a1_vec[l] = a;
            if (cumul >= rando) {
              break;
            }
          }
        }

        if(I[I_dx(l, sim_missing[i]-1, 1, 2, N)] == 0) {
          a2_vec[l] = -1;
        } else {
          rando = runif(1)[0] * S2;
          cumul = 0.0;
          for(a = 0; a < A[l]; a++) {
            cumul += AC[D_dx(l, sim_colls[i]-1, a, L, C, A, CA)] - (a1_vec[l] == a);
            a2_vec[l] = a;
            if (cumul >= rando) {
              break;
            }
          }
        }
      }
    }
    for(c = 0; c < C; c++) {
      sum = 0.0;
      LOO = c == (sim_colls[i]-1);
      for(l = 0; l < L; l++) {
        GPROB_DIP_FROM_SIM(a1_vec[l], a2_vec[l], l, c, gp);
        sum += log(gp);
      }
      out(c, i) = sum;
    }
  }
  return(out);
}
