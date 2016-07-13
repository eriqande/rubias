#include <Rcpp.h>
#include "macros.h"
#include "utilities.h"
using namespace Rcpp;

//' Simulate genotype log-likelihoods from a population by gene copy
//'
//' @param par_list genetic data converted to the param_list format by tcf2param_list
//'
//' @param sim_colls a vector of indices for the collections desired for simulation;
//' each element of the list corresponds to an individual
//'
//'
//' @return \code{gprob_sim} returns a matrix of the summed log-likelihoods
//' for all loci of a simulated population mixture; columns represent individuals,
//' with each row containing their log-likelihood of belonging to the collection
//' of the same index, given the selection of two independent gene copies from the desired
//' collection of origin's reference allele frequencies
//'
//' @examples
//' example(tcf2param_list)
//' sim_colls <- sample(ale_par_list$C, 1070, replace = T)
//' ale_sim_gprobs_gc <- gprob_sim_gc(ale_par_list, sim_colls)
//' ale_sim_gprobs_ind <- gprob_sim_ind(ale_par_list, sim_colls)
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
//' @param par_list genetic data converted to the param_list format by tcf2param_list
//'
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
//' @examples
//' example(tcf2param_list)
//' sim_colls <- sample(ale_par_list$C, 1070, replace = T)
//' ale_sim_gprobs <- gprob_sim(ale_par_list, sim_colls)
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
  double rando, sum, gp;
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

//' Simulate genotype log-likelihoods from a population by gene copy
//'
//' @param par_list genetic data converted to the param_list format by tcf2param_list
//'
//' @param sim_colls a vector of indices for the collections desired for simulation;
//' each element of the list corresponds to an individual
//'
//'
//' @return \code{gprob_sim} returns a matrix of the summed log-likelihoods
//' for all loci of a simulated population mixture; columns represent individuals,
//' with each row containing their log-likelihood of belonging to the collection
//' of the same index, given the selection of two independent gene copies from the desired
//' collection of origin's reference allele frequencies
//'
//' @examples
//' example(tcf2param_list)
//' sim_colls <- sample(ale_par_list$C, 1070, replace = T)
//' ale_sim_gprobs_gc <- gprob_sim_gc(ale_par_list, sim_colls)
//' ale_sim_gprobs_ind <- gprob_sim_ind(ale_par_list, sim_colls)
//'
//' @export
// [[Rcpp::export]]
List gprob_sim_gc_RU(List par_list, IntegerVector sim_colls, IntegerVector ru_colls) {
  int i, ru, c, l , a, S, S2;
  int N = sim_colls.size();
  int C = as<int>(par_list["C"]);
  int K = as<int>(par_list["K"]);
  int L = as<int>(par_list["L"]);
  int LOO;
  IntegerVector a1_vec(L);
  IntegerVector a2_vec(L);
  IntegerVector AC = as<IntegerVector>(par_list["AC"]);
  IntegerVector AC_R = as<IntegerVector>(par_list["AC_R"]);
  IntegerVector sum_AC = as<IntegerVector>(par_list["sum_AC"]);
  IntegerVector sum_AC_R = as<IntegerVector>(par_list["sum_AC_R"]);
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector DP_R = as<NumericVector>(par_list["DP_R"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  NumericVector sum_DP_R = as<NumericVector>(par_list["sum_DP_R"]);
  NumericVector RU_starts = as<NumericVector>(par_list["RU_starts"]);
  NumericVector RU_vec = as<NumericVector>(par_list["RU_vec"]);
  double cumul, rando, sum, gp;
  NumericMatrix out(C,N);
  NumericMatrix out_r(K,N);


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
    for(ru = 0; ru < K; ru++) {
      sum = 0.0;
      LOO = ru == (ru_colls[i]-1);
      for(l = 0; l < L; l++){
        GPROB_DIP_FROM_SIM_R(a1_vec[l], a2_vec[l], l, ru, gp);
        sum += log(gp);
      }
      out_r(ru, i) = sum;
    }
  }
  return(List::create(out, out_r));
}

//' Simulate genotype log-likelihoods from a population by individual for Ben's experiment
//'
//' @param par_list genetic data converted to the param_list format by tcf2param_list
//'
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
//' @examples
//' example(tcf2param_list)
//' sim_colls <- sample(ale_par_list$C, 1070, replace = T)
//' ale_sim_gprobs <- gprob_sim(ale_par_list, sim_colls)
//'
//' @export
// [[Rcpp::export]]
List gprob_sim_ind_RU(List par_list, IntegerVector sim_colls) {
  int i, c, l, k, n, ru, ind, count;
  int Ind = sim_colls.size();
  int N = as<int>(par_list["N"]);
  int K = as<int>(par_list["K"]);
  int C = as<int>(par_list["C"]);
  int L = as<int>(par_list["L"]);
  int LOO;
  int sim_ru;
  IntegerVector a1_vec(L);
  IntegerVector a2_vec(L);
  IntegerVector coll_N = as<IntegerVector>(par_list["coll_N"]);
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  IntegerVector I_coll = as<IntegerVector>(par_list["coll"]);
  IntegerVector A = as<IntegerVector>(par_list["A"]);
  IntegerVector CA = as<IntegerVector>(par_list["CA"]);
  NumericVector DP = as<NumericVector>(par_list["DP"]);
  NumericVector DP_R = as<NumericVector>(par_list["DP_R"]);
  NumericVector sum_DP = as<NumericVector>(par_list["sum_DP"]);
  NumericVector sum_DP_R = as<NumericVector>(par_list["sum_DP_R"]);
  NumericVector RU_starts = as<NumericVector>(par_list["RU_starts"]);
  NumericVector RU_vec = as<NumericVector>(par_list["RU_vec"]);
  double rando, sum, gp;
  NumericMatrix out(C,Ind);
  NumericMatrix out_r(K,Ind);


  for(i = 0; i < Ind; i++) {
    // rebuild an RUs list from sim_colls
    for(ru = 0; ru < K; ru++) {
      for(c = RU_starts[ru]; c < RU_starts[ru + 1]; c++) {
        if(RU_vec[c] == sim_colls[i]) {
          break;
        }
      }
      sim_ru = ru + 1;
    }


    n = coll_N[sim_colls[i] - 1];
    ind = randint(n);
    count = 0;
    for(k = 0; k < N; k++) {           // Not a terribly efficient way to get the (ind)th individual of a collection, but best way without knowing I is organized by collection, and without a new indexing macro
      if(I_coll[k] == sim_colls[i]) {
        count += 1;
      }
      if(ind == count) {
        break;
      }
    }
    for(l = 0; l < L; l++) {
      a1_vec[l] = I[I_dx(l, k, 0, 2, N)] - 1;
      a2_vec[l] = I[I_dx(l, k, 1, 2, N)] - 1;
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
    for(ru = 0; ru < K; ru++) {
      sum = 0.0;
      LOO = c == (sim_ru - 1);
      for(l = 0; l < L; l ++) {
        GPROB_DIP_FROM_SIM_R(a1_vec[l], a2_vec[l], l, ru, gp);
      }
      out_r(ru,i) = sum;
    }
  }
  return(List::create(out, out_r));
}

//' Simulate genotypes by gene copy, with missing data from chosen individuals
//'
//' @param par_list list_diploid_params output
//' @param sim_colls a vector length I of collections from which to sample the genotypes
//' for individual i
//' @param sim_missing a vector of length I of indices for individuals in I_list
//' whose missing data should be copied for individual i
//'
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
          printf("l = %d, indiv = %d \n", l, sim_missing[i]);
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
