#include <RcppArmadillo.h>
#include <splines.h>
#include <c_optim.h>

namespace rstpm2 {

  class gsm_functor {
  public:
    double inflate, target;
    arma::vec gamma;
    double eps;
    ns ns1;
    double link(double S) { return std::log(-std::log(S)); }
    double linkinv(double eta) { return std::exp(-std::exp(eta)); }
    gsm_functor(arma::vec knots, arma::vec Boundary_knots, int intercept, arma::vec gamma,
		arma::mat q_const,
		double inflate = 100.0, double target = 0.0, double eps = 1.0e-12) :
      inflate(inflate), target(target), gamma(gamma), eps(eps) {
      // ns(vec boundary_knots, vec interior_knots, mat _q_matrix, int intercept=0, int cure=0) :
      ns1 = ns(Boundary_knots, knots, q_const, intercept, 0);
    }
    double eta(double t, double etap) {
      return arma::sum(ns1.eval(std::log(t),0) % gamma) + etap;
    }
    double operator()(double t) {
      return eta(t,0.0) - target;
    }
    double rand(double etap, double tentry, double maxt) {
      double u = R::runif(0.0,1.0);
      target = (tentry==0.0 ? link(u) : link(u*linkinv(eta(tentry,etap)))) - etap;
      return std::get<0>(R_zeroin2_functor_ptr<gsm_functor>(std::max(tentry/inflate,eps),
							    maxt*inflate, this, 1.0e-10));
    }
  };

  RcppExport SEXP call_rgsm(SEXP n_, SEXP tentry_, SEXP maxt_, SEXP etap_, SEXP knots_,
			    SEXP Boundary_knots_, SEXP intercept_, SEXP gamma_, SEXP q_const_,
			    SEXP inflate_ /*= 100.0*/,
			    SEXP eps_ /*= 1.0e-12*/) {
    Rcpp::RNGScope rngScope;
    int n = Rcpp::as<int>(n_);
    double tentry = Rcpp::as<double>(tentry_); // arma::vec?
    double maxt = Rcpp::as<double>(maxt_);
    double etap = Rcpp::as<double>(etap_);
    arma::vec knots = Rcpp::as<arma::vec>(knots_);
    arma::vec Boundary_knots = Rcpp::as<arma::vec>(Boundary_knots_);
    int intercept = Rcpp::as<int>(intercept_);
    arma::vec gamma = Rcpp::as<arma::vec>(gamma_);
    arma::mat q_const = Rcpp::as<arma::mat>(q_const_);
    double inflate = Rcpp::as<double>(inflate_);
    double eps = Rcpp::as<double>(eps_);
    //
    std::vector<double> y(n);
    gsm_functor gsm(knots, Boundary_knots, intercept, gamma, q_const, inflate, eps);
    for (int i=0; i<n; i++)
      y[i] = gsm.rand(etap, tentry, maxt);
    return Rcpp::wrap(y);
  }
  
} // namespace rstpm2
