#include <RcppArmadillo.h>
#include <c_optim.h>
#include <splines.h>

namespace rstpm2 {
  using namespace arma;
  using namespace Rcpp;

  class aft {
  public:
    List args;
    vec init;
    mat X, X0;
    mat XD, XD0;
    vec event;
    vec time, time0;
    vec boundaryKnots;
    vec interiorKnots;
    mat q_matrix;
    int cure;
    ns s;
    bool delayed;
    double kappa; // scale for the quadratic penalty for monotone splines
    aft(SEXP list);
    mat rmult(mat m, vec v);
    mat rmult(mat m, uvec v);
    double objective(vec betafull);
    vec gradientPenalty(mat Q, vec beta);
    vec gradient(vec betafull);
    double objective(NumericVector betafull);
    NumericVector gradient(NumericVector betafull);
    vec survival(vec time, mat X);
    vec haz(vec time, mat X, mat XD);
    mat gradh(vec time, mat X, mat XD);
  };
  

} // namespace rstpm2

