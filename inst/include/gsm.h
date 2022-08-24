#ifndef RSTPM2_GSM_H
#define RSTPM2_GSM_H

#include <RcppArmadillo.h>

namespace rstpm2 {

  enum link_types {PH};
  
  class gsm {
  public:
    link_types link_type;
    double tmin, tmax, target;
    arma::vec gamma, etap;
    bool log_time;
    ns ns1;
    double link(double S);
    double linkinv(double eta);
    gsm(); // default constructor
    gsm(std::string link_name, arma::vec knots, arma::vec Boundary_knots, int intercept,
	arma::vec gamma, arma::vec etap,
	arma::mat q_const, int cure,
	double tmin = 1.0e-4, double tmax = 1.0e4, double inflate = 100.0, bool log_time = true);
    double eta(double y, double etap);
    double operator()(double y);
    double rand(double etap, double tentry=0.0);
  };

  gsm read_gsm(SEXP args);
  
} // namespace rstpm2

#endif /* gsm.h */
