#ifndef C_OPTIM_H
#define C_OPTIM_H

#include <Rcpp.h>

namespace rstpm2 {

  typedef double optimfn(int, double *, void *);
  typedef void optimgr(int, double *, double *, void *);

  double min(double a, double b);
  double max(double a, double b);
  double bound(double x, double lower, double upper);

  class NelderMead {
  public:
    NelderMead(int trace = 0, int maxit = 500, 
	       double abstol = - INFINITY,
	       double reltol = 1.0e-8, 
	       double alpha = 1.0, double beta = 0.5, double gamma = 2.0);
    void optim(optimfn fn, Rcpp::NumericVector init, void * ex);
    int n, trace, maxit, fail, fncount;
    double abstol, reltol, alpha, beta, gamma, Fmin;
    Rcpp::NumericVector coef;
  };

  class BFGS {
  public:
    BFGS(int trace = 0, int maxit = 100, 
	 double abstol = - INFINITY,
	 double reltol = 1.0e-8, int report = 10);
    void optim(optimfn fn, optimgr gr, Rcpp::NumericVector init, void * ex,
	       double eps = 1.0e-8);
    double calc_objective(optimfn fn, Rcpp::NumericVector coef, void * ex);
    double calc_objective(optimfn fn, void * ex);
    Rcpp::NumericMatrix calc_hessian(optimgr gr, void * ex, double eps = 1.0e-8);
    int n, trace, maxit, report, fncount, grcount, fail;
    double abstol, reltol, Fmin;
    Rcpp::NumericVector coef;
    Rcpp::NumericMatrix hessian;
  };

  typedef double (*Brent_fminfn)(double, void *);

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		  void *info, double tol);

} // anonymous rstpm2

#endif /* c_optim_h */
