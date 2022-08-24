#include <RcppArmadillo.h>
#include <splines.h>
#include <c_optim.h>
#include <gsm.h>

namespace rstpm2 {

  double gsm::link(double S) {
    return link_type==PH ? std::log(-std::log(S)) : -100.0;
  }
  double gsm::linkinv(double eta) {
    return link_type==PH ? std::exp(-std::exp(eta)) : 1.0e-10;
  }
  gsm::gsm() {}
  gsm::gsm(std::string link_name,
	   arma::vec knots, arma::vec Boundary_knots, int intercept, arma::vec gamma,
	   arma::vec etap,
	   arma::mat q_const, int cure,
	   double tmin, double tmax, double inflate, bool log_time) :
    tmin(tmin/inflate), tmax(tmax*inflate), gamma(gamma), etap(etap), log_time(log_time) {
    // ns(vec boundary_knots, vec interior_knots, mat _q_matrix, int intercept=0, int cure=0) :
    ns1 = ns(Boundary_knots, knots, q_const, intercept, cure);
    target = 0.0;
    if (link_name == "PH") link_type = PH;
  }
  double gsm::eta(double y, double etap) {
    return arma::sum(ns1.eval(y,0) % gamma) + etap;
  }
  double gsm::operator()(double y) {
    return eta(y,0.0) - target;
  }
  double gsm::rand(double etap, double tentry) {
    using std::log;
    double u = R::runif(0.0,1.0);
    double ymin = tentry == 0.0 ? (log_time ? log(tmin) : tmin) : (log_time ? log(tentry) : tentry);
    double ymax = log_time ? log(tmax) : tmax;
    target = (tentry==0.0 ? link(u) : link(u*linkinv(eta(ymin, etap)))) - etap;
    double root = std::get<0>(R_zeroin2_functor_ptr<gsm>(ymin, ymax, this, 1.0e-8, 100));
    return log_time ? std::exp(root) : root;
  }

  gsm read_gsm(SEXP list_) {
    try {
      using Rcpp::as;
      Rcpp::List list = as<Rcpp::List>(list_);
      std::string link_name = as<std::string>(list("link_name"));
      double tmin = as<double>(list("tmin"));
      double tmax = as<double>(list("tmax"));
      arma::vec knots = as<arma::vec>(list("knots"));
      arma::vec Boundary_knots = as<arma::vec>(list("Boundary_knots"));
      int intercept = as<int>(list("intercept"));
      arma::vec gamma = as<arma::vec>(list("gamma"));
      arma::vec etap = as<arma::vec>(list("etap"));
      arma::mat q_const = as<arma::mat>(list("q_const"));
      int cure = as<int>(list("cure"));
      double inflate = as<double>(list("inflate"));
      bool log_time = as<bool>(list("log_time"));
      return gsm(link_name, knots, Boundary_knots, intercept, gamma,
		 etap,
		 q_const, cure, tmin, tmax, inflate, log_time);
    } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
    } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
    }
    return gsm();
  } 

  RcppExport SEXP test_read_gsm(SEXP args) {
    Rcpp::RNGScope rngScope;
    gsm gsm1 = read_gsm(args);
    return Rcpp::wrap(gsm1.rand(gsm1.etap(0)));
  }
  
} // namespace rstpm2
