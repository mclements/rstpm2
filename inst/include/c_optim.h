#ifndef C_OPTIM_H
#define C_OPTIM_H

#include <Rcpp.h>

namespace rstpm2 {

  typedef double optimfn(int, double *, void *);
  typedef void optimgr(int, double *, double *, void *);

  /* type of pointer to the target and gradient functions for Nlm */
  typedef void (*fcn_p)(int, double *, double *, void *);

  /* type of pointer to the hessian functions for Nlm */
  typedef void (*d2fcn_p)(int, int, double *, double *, void *);

  double min(double a, double b);
  double max(double a, double b);
  double bound(double x, double lower, double upper);

  /**
     Adapt a function object (functor) for NelderMead and BFGS
  **/
  template<class T>
    double adapt_functor(int n, double * beta, void * par) {
    T * model = (T *) par;
    Rcpp::NumericVector x(beta,beta+n);
    return model->operator()(x);
  }
  /**
     Adapt an objective function for NelderMead and BFGS
  **/
  template<class T>
    double adapt_objective(int n, double * beta, void * par) {
    T * model = (T *) par;
    Rcpp::NumericVector x(beta,beta+n);
    return model->objective(x);
  }
  /**
     Adapt a gradient function for BFGS
  **/
  template<class T>
    void adapt_gradient(int n, double * beta, double * grad, void * par) {
    T * model = (T *) par;
    Rcpp::NumericVector x(beta,beta+n);
    Rcpp::NumericVector vgrad = model->gradient(x);
    for (int i=0; i<n; ++i) grad[i] = vgrad[i];
  }

  class NelderMead {
  public:
    NelderMead(int trace = 0, int maxit = 500, 
	       double abstol = - INFINITY,
	       double reltol = 1.0e-8, 
	       double alpha = 1.0, double beta = 0.5, double gamma = 2.0, 
	       double epshess = 6.055454e-06, bool hessianp = true);
    virtual void optim(optimfn fn, Rcpp::NumericVector init, void * ex);
    template<class T>
      void optim(Rcpp::NumericVector init, T object) {
      optim(&adapt_objective<T>,init,(void *) &object);
    }
    virtual Rcpp::NumericMatrix calc_hessian(optimfn fn, void * ex);
    int n, trace, maxit, fail, fncount;
    double abstol, reltol, alpha, beta, gamma, Fmin, epshess;
    bool hessianp;
    Rcpp::NumericVector coef;
    Rcpp::NumericMatrix hessian;
  };

  class BFGS {
  public:
    BFGS(int trace = 0, int maxit = 100, 
	 double abstol = - INFINITY,
	 double reltol = 1.0e-8, int report = 10, double epshess = 1.0e-8, bool hessianp = true);
    virtual void optim(optimfn fn, optimgr gr, Rcpp::NumericVector init, void * ex);
    virtual void optim(int n, optimfn fn, optimgr gr, double * init, void * ex);
    virtual double calc_objective(optimfn fn, Rcpp::NumericVector coef, void * ex);
    virtual double calc_objective(optimfn fn, void * ex);
    virtual Rcpp::NumericMatrix calc_hessian(optimgr gr, void * ex);
    template<class T>
      void optim(Rcpp::NumericVector init, T object) {
      optim(&adapt_objective<T>,&adapt_gradient<T>,init,(void *) &object);
    }
    int n, trace, maxit, report, fncount, grcount, fail;
    double abstol, reltol, Fmin, epshess;
    bool hessianp;
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
  double epshess = 6.055454e-06,
	int itrmcd = 0, 
	int itncnt = 0,
	bool hessianp = true);
    void optim(fcn_p fcn, fcn_p d1fcn, Rcpp::NumericVector init, void * state);
    void optim(fcn_p fcn, Rcpp::NumericVector init, void * state);
    double calc_objective(fcn_p fn, Rcpp::NumericVector coef, void * ex);
    double calc_objective(fcn_p fn, void * ex);
    Rcpp::NumericMatrix calc_hessian(fcn_p gr, void * ex);
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
    double epshess;
    int itrmcd;
    int itncnt;
    bool hessianp;
    Rcpp::NumericVector coef;
    Rcpp::NumericMatrix hessian;
  };


  typedef double (*Brent_fminfn)(double, void *);

  double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		    void *info, double tol);

  double R_zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit);				/* Max # of iterations */

/** 
    Adapt a function object (functor) to work with R_zeroin2()
**/
template<class Functor>
double R_zeroin2_adaptor(double x, void * par) {
  Functor * functor = (Functor *) par;
  return functor->operator()(x);
}

/** 
    Use R_zeroin2 with a function object (functor)
    @return tuple<double,double,int> with (root, Tol, Maxit)
**/
template<class Functor>
std::tuple<double,double,int>
R_zeroin2_functor(double a, double b, Functor functor, double eps = 1.0e-8) {
  double Tol = eps;
  int Maxit = 100;
  double root = R_zeroin2(a,b,functor(a),functor(b),&R_zeroin2_adaptor<Functor>,(void *) &functor,
			  &Tol, &Maxit);
  return std::make_tuple(root, Tol, Maxit);
}

/** 
    Use R_zeroin2 with a function object pointer (functor)
    @return tuple<double,double,int> with (root, Tol, Maxit)
**/
template<class Functor>
std::tuple<double,double,int>
R_zeroin2_functor_ptr(double a, double b, Functor *functor, double tol = 1.0e-8, int maxit = 100) {
  double Tol = tol;
  int Maxit = maxit;
  double root = R_zeroin2(a,b,(*functor)(a),(*functor)(b),
			  &R_zeroin2_adaptor<Functor>,(void *) functor,
			  &Tol, &Maxit);
  return std::make_tuple(root, Tol, Maxit);
}

  
  
  /** 
      Adapt a function object (functor) to work with Brent_fmin()
  **/
  template<class T, class X>
    double Brent_fmin_functor(X x, void * par) {
    T * model = (T *) par;
    return model->operator()(x);
  }

  /** 
      Use Brent_fmin with a function object (functor)
  **/
  template<class T>
    double BrentFmin(double a, double b, T obj, double eps = 1.0e-8) {
    return Brent_fmin(a,b,&Brent_fmin_functor<T,double>,(void *) &obj,eps);
  }

  Rcpp::NumericMatrix qr_q(const Rcpp::NumericMatrix& X, double tol = 1E-12); 

  
} // anonymous rstpm2

#endif /* c_optim_h */
