#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(rTRNG, RcppParallel)]]

#include <RcppParallel.h>
using namespace RcppParallel;
#include <trng/yarn2.hpp>
#include <trng/normal_dist.hpp>
using namespace trng;


typedef yarn2 rngKind;


struct KernelWorker : public Worker {

  // Input data
  const int J, K;
  const RVector<int> j, mk, agg;
  const RVector<double> V0, R, beta, PD;
  const RMatrix<double> Z, r;
  const std::vector<rngKind> rngSplit; // pre-split sub-sequences
  // Output object to write to
  RMatrix<double> L_agg;

  // Initialize from input and output objects (automatic wrap of Rcpp objects
  // as RMatrix/RVector)
  KernelWorker(const NumericVector V0, const NumericVector R,
               const NumericVector beta,  const NumericVector PD,
               const IntegerVector j, const int J,
               const NumericMatrix Z, const NumericMatrix r,
               const int K, const IntegerVector mk,
               const IntegerVector agg,
               const std::vector<rngKind> rngSplit,
               NumericMatrix L_agg)
    : J(J), K(K), j(j), mk(mk), agg(agg), V0(V0), R(R), beta(beta), PD(PD),
      Z(Z), r(r), rngSplit(rngSplit), L_agg(L_agg)
  {}

  // operator processing an exclusive range of indices
  void operator()(std::size_t begin, std::size_t end) {

    rngKind rngj;
    normal_dist<> normal(0.0, 1.0); // TRNG normal distribution

    unsigned int mk_begin = (unsigned int)begin;
    unsigned int mk_end = (unsigned int)end;
    int thismk, prevmk;
    unsigned int m;

    int thisj;
    double sigmaj, thetaj, Yj, Lj;

    // loop over the set of relevant counterparties j
    for (unsigned int j_idx = 0; j_idx < j.length(); j_idx++) {

      thisj = j[j_idx] - 1; // current counterparty (0-based indexing)
      sigmaj = sqrt(1 - beta[j_idx]*beta[j_idx]);
      thetaj = normal.icdf(PD[j_idx]);
      rngj = rngSplit[j_idx]; // j-th subsequence

      prevmk = -1;
      // loop over the relevant combined simulation indices
      for (unsigned int mk_idx = mk_begin; mk_idx < mk_end; mk_idx++) {
        thismk = mk[mk_idx] - 1; // current combined simulation
        rngj.jump(thismk - prevmk - 1); // jump to thismk from prevmk
        m = thismk/K; // market scenario simulation index

        // credit environment combined return
        Yj = beta[j_idx]*Z(m, thisj) + sigmaj*normal(rngj);
        // default indicator and corresponding loss
        Lj = V0[j_idx] - ( Yj < thetaj ?
                             R[j_idx] :
                             V0[j_idx] * (1 + r(m, thisj)) );
        L_agg(mk_idx, agg[j_idx] - 1) += Lj;

        prevmk = thismk;
      }

    }

  }

};

// [[Rcpp::export]]
void simulationKernel_C(const DataFrame pf,
                        const NumericMatrix Z, const NumericMatrix r,
                        const int J,
                        const int K, const IntegerVector mk,
                        const IntegerVector agg,
                        const unsigned long seed,
                        NumericMatrix L_agg) {

  // counterparty indices in the (sub)portfolio
  IntegerVector j = pf["j"];
  // RNG for the full sequence
  rngKind rngFull(seed);
  // RNGs for the LeapFrog subsequences of individual counterparties
  std::vector<rngKind> rngSplit(j.length());
  for (unsigned int j_idx = 0; j_idx < j.length(); j_idx++) {
    int thisj = j[j_idx] - 1; // 0-based C++ indexing
    // split the full sequence based on the full set of counterparties (J)
    rngSplit[j_idx] = rngFull;
    rngSplit[j_idx].split(J, thisj);
  }

  // worker for the parallel simulation
  KernelWorker w(pf["V0"], pf["R"], pf["beta"], pf["PD"],
                 j, J,
                 Z, r,
                 K, mk,
                 agg,
                 rngSplit,
                 L_agg);

  // parallel execution over the set of relevant simulation indices mk
  int grainSize = 100;
  parallelFor(0, mk.length(), w, grainSize);

}
