#include <Rcpp.h>
#include "macros.h"
using namespace Rcpp;

//' Find all pairs that have close matches
//'
//'
//' @keywords internal
//'
//' @param par_list genetic data converted to the param_list format by \code{tcf2param_list}
//'
//' @return Gotta say more
//' @examples
//' # gotta do something here too
//' @export
// [[Rcpp::export]]
List rcpp_close_matchers(List par_list, double non_miss_fract, double match_fract) {
  int i, j, l, nm, ma;
  int q0, q1, r0, r1; // to store first and second gene copy from first indiv (q) and second (r)
  IntegerVector I = as<IntegerVector>(par_list["I"]);
  int N = as<int>(par_list["N"]);
  int L = as<int>(par_list["L"]);
  std::vector<int> nm_vec; // num-non-missing vector
  std::vector<int> ma_vec; // num-matching vector
  std::vector<int> i_vec; // index of first fish vector
  std::vector<int> j_vec; // index of first fish vector

  for(i = 0; i < N; i++) { // cycle over first individual
    for(j = i + 1; j < N; j++) {  // cycle over second individual
      nm = 0;
      ma = 0;
      for(l = 0; l < L; l++) {  // cycle over loci.  Note the assumption that all loci are diploid
        q0 = I[I_dx(l, i, 0, 2, N)] - 1;  // all these are base 0.  -1 means missing
        q1 = I[I_dx(l, i, 1, 2, N)] - 1;
        r0 = I[I_dx(l, j, 0, 2, N)] - 1;
        r1 = I[I_dx(l, j, 1, 2, N)] - 1;

        if(q0 > -1 && q1 > -1 && r0 > -1 && r1 > -1) {
          nm++;  // increment number of non-missing loci
          if( (q0 == r0 && q1 == r1) || (q0 == r1 && q1 == r0)) {  // if the genotypes match
            ma++;
          }
        }

      }

      if( (double)nm / L >= non_miss_fract   &&    (double)ma/nm >= match_fract ) {
        i_vec.push_back(i + 1);
        j_vec.push_back(j + 1);
        nm_vec.push_back(nm);
        ma_vec.push_back(ma);
      }
    }
  }

  return List::create(
    _["indx1"] = i_vec,
    _["indx2"] = j_vec,
    _["num_non_miss"] = nm_vec,
    _["num_match"] = ma_vec
  );
}
