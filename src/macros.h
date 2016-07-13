// Indexing macros
// i = individual index (from 0)
// l = locus index
// g = gene copy index
// ru = reporting unit index
// c = collection index
// a = allele index
// N = Total # of individuals
// L = Total # of loci
// C = Total # of collections
// P = Ploidy of individuals
// A = Vector of number of alleles at each locus
// CA = Vector of cumulative number of alleles at all previous loci
// AC = Vector of allele counts (Locus -> Population -> allele)
// sum_AC = Vector of sums of allele counts at each locus (Loci -> Population)
// DP = Vector of Dirichlet parameters (Locus -> Population -> allele)
// I = Vector of Individual genotypes (Locus -> Individual -> allele)
// LOO = Leave-one-out cross validation? 0 or 1
// sum_DP = vector of sums of Dirichlet parameters (Loci -> Population)
// a1 = allelic type 1 from simulation (base 0)
// a2 = allelic type 2 from simulation (base 0)
#define D_dx(l, c, a, L, C, A, CA) (C) * (CA[l]) + (A[l]) * (c) + (a)
#define I_dx(l, i, g, P, N) (P) * (N) * (l) + (P) * (i) + (g)
#define SD_dx(l, c, C) (C) * (l) + (c)
#define GPROB_DIP(i, l, c, result) {                                                                  \
int a1 = I[I_dx(l, i, 0, 2, N)] - 1;                                                                  \
int a2 = I[I_dx(l, i, 1, 2, N)] - 1;                                                                  \
if(a1 < 0 || a2 < 0) {result = 1.0;} else {                                                           \
  double y1 = DP[D_dx(l, c, a1, L, C, A, CA)];                                                        \
  double y2 = DP[D_dx(l, c, a2, L, C, A, CA)];                                                        \
  double n = sum_DP[SD_dx(l, c, C)];                                                                  \
  double nsub = LOO * (1 + (a1 == a2));                                                               \
  double dsub = LOO * 2;                                                                              \
  result = (y1 - nsub) * (y2 - nsub + (a1 == a2)) * (1 + (a1 != a2)) / ((n - dsub) * (n + 1 - dsub)); \
}                                                                                                     \
}
#define GPROB_DIP_FROM_SIM(a1, a2, l, c, result) {                                                                  \
if(a1 < 0 || a2 < 0) {result = 1.0;} else {                                                                       \
  double y1 = DP[D_dx(l, c, a1, L, C, A, CA)];                                                                      \
  double y2 = DP[D_dx(l, c, a2, L, C, A, CA)];                                                                      \
  double n = sum_DP[SD_dx(l, c, C)];                                                                                \
  double nsub = LOO * (1 + (a1 == a2));                                                                             \
  double dsub = LOO * 2;                                                                                            \
  result = (y1 - nsub) * (y2 - nsub + (a1 == a2)) * (1 + (a1 != a2)) / ((n - dsub) * (n + 1 - dsub));               \
  if(result <= 0 || std::isnan(result)) printf("y1 = %f y2 = %f result = %f a1 = %d a2 = %d i = %d  l = %d  c = %d nsub = %f dsub = %f rando = %f \n", y1, y2, result, a1, a2, i, l, c, nsub, dsub, rando);\
}                                                                                                                   \
}                                                              \

#define GPROB_DIP_FROM_SIM_R(a1, a2, l, ru, result) {                                                                                                                                                         \
if(a1 < 0 || a2 < 0) {result = 1.0;} else {                                                                                                                                                                \
  double y1 = DP_R[D_dx(l, ru, a1, L, K, A, CA)];                                                                                                                                                             \
  double y2 = DP_R[D_dx(l, ru, a2, L, K, A, CA)];                                                                                                                                                             \
  double n = sum_DP_R[SD_dx(l, ru, K)];                                                                                                                                                                       \
  double nsub = LOO * (1 + (a1 == a2));                                                                                                                                                                    \
  double dsub = LOO * 2;                                                                                                                                                                                   \
  result = (y1 - nsub) * (y2 - nsub + (a1 == a2)) * (1 + (a1 != a2)) / ((n - dsub) * (n + 1 - dsub));                                                                                                      \
  if(result <= 0 || std::isnan(result)) printf("y1 = %f y2 = %f result = %f a1 = %d a2 = %d i = %d  l = %d  c = %d nsub = %f dsub = %f rando = %f \n", y1, y2, result, a1, a2, i, l, c, nsub, dsub, rando);\
}                                                                                                                                                                                                          \
}
