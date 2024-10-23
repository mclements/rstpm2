#include <RcppArmadillo.h>
#include <c_optim.h>
#include <splines.h>

namespace rstpm2 {
  using namespace arma;
  using namespace Rcpp;

  class aft : public ConstrBFGSx {
  public:
    List args;
    vec init;
    std::vector<mat> X_vector, X_vector0;
    mat Xc;
    mat X, X0;
    mat XD, XD0;
    vec event;
    vec time, time0;
    vec gnodes, gweights, bhazard;
    vec boundaryKnots;
    vec interiorKnots;
    mat q_matrix;
    int cure;
    int mixture, tvc_integrated;
    ns s;
    double kappa; // scale for the quadratic penalty for monotone splines
    double eps1;  // minimum hazard
    uvec delayed;
    bool add_penalties;
    aft(SEXP list);
    mat rmult(mat m, vec v);
    mat rmult(mat m, uvec v);
    double objective(vec betafull);
    double objective_integrated(vec betafull);
    double objective_cumulative(vec betafull);
    vec gradientPenalty(mat Q, vec beta);
    vec gradient(vec betafull);
    vec gradient_integrated(vec betafull);
    vec gradient_cumulative(vec betafull);
    double objective(NumericVector betafull);
    NumericVector gradient(NumericVector betafull);
    // vec survival(vec time, mat X);
    // vec haz(vec time, mat X, mat XD);
    // mat gradh(vec time, mat X, mat XD);
  };
  

} // namespace rstpm2

