#ifndef C_OPTIM_H
#define C_OPTIM_H

#include <Rcpp.h>

namespace rstpm2 {

  typedef double optimfn(int, double *, void *);
  typedef void optimgr(int, double *, double *, void *);

  /* type of pointer to the target and gradient functions  for Nlm */
  typedef void (*fcn_p)(int, double *, double *, void *);

  /* type of pointer to the hessian functions */
  typedef void (*d2fcn_p)(int, int, double *, double *, void *);

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
    virtual void optim(optimfn fn, optimgr gr, Rcpp::NumericVector init, void * ex,
	       double eps = 1.0e-8);
    virtual double calc_objective(optimfn fn, Rcpp::NumericVector coef, void * ex);
    virtual double calc_objective(optimfn fn, void * ex);
    virtual Rcpp::NumericMatrix calc_hessian(optimgr gr, void * ex, double eps = 1.0e-8);
    int n, trace, maxit, report, fncount, grcount, fail;
    double abstol, reltol, Fmin;
    Rcpp::NumericVector coef;
    Rcpp::NumericMatrix hessian;
  };

  class Nlm {
  public:
    Nlm(double fscale = 1.0,    // nlm()
	int method = 2,         // cf. nlm: method=1
	int iexp = 1,           // nlm()
	int msg = 9,            // nlm()
	int ndigit = 12,        // nlm()
	int itnlim = 50,        // nlm()
	int iagflg = 1,         // nlm()
	int iahflg = 0,         // nlm()
	double dlt = 1.0,       // nlm
	double gradtl = 1.0e-6, // nlm()
	double stepmx = 0.0,    // set to -1.0 to get nlm()'s behaviour
	double steptl = 1.0e-6,  // nlm()
	int itrmcd = 0, 
	int itncnt = 0
	);
    void optim(fcn_p fcn, fcn_p d1fcn, Rcpp::NumericVector init, void * state,
	       bool hessianp = true); // assumes iahflg=0
    double calc_objective(fcn_p fn, Rcpp::NumericVector coef, void * ex);
    double calc_objective(fcn_p fn, void * ex);
    Rcpp::NumericMatrix calc_hessian(fcn_p gr, void * ex, double eps = 1.0e-8);
    void set_print_level(int);
    double fscale;
    int method;
    int iexp;
    int msg;
    int ndigit;
    int itnlim;
    int iagflg;
    int iahflg;
    double dlt;
    double gradtl;
    double stepmx;
    double steptl;
    int itrmcd;
    int itncnt;
    Rcpp::NumericVector coef;
    Rcpp::NumericMatrix hessian;
  };


  typedef double (*Brent_fminfn)(double, void *);

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		  void *info, double tol);

} // anonymous rstpm2

#endif /* c_optim_h */
