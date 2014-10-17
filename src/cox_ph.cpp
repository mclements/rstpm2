#include <RcppArmadillo.h>
#include "c_optim.h"

namespace rstpm2 {

  using namespace Rcpp;
  using namespace arma;

  // assume right censored values in ascending order of time and no ties

  RcppExport SEXP test_cox_tvc2(SEXP args) {
  
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    vec x      = as<vec>(largs["x"]); // one covariate, length n
    vec beta   = as<vec>(largs["beta"]); // length 2
    int n      = time.size();
  
    double llike = 0.0;
    vec eta;
  
    for (int i=0; i<n; ++i) {
      if (event(i) == 1) {
	eta = beta(0)*x(span(i,n-1)) + beta(1)*log(time(i))*x(span(i,n-1));
	llike += eta(0) - log(sum(exp(eta)));
      }
    }
    return wrap(llike);
  }

  RcppExport SEXP test_cox_tvc2_grad(SEXP args) {
  
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    vec x      = as<vec>(largs["x"]); // one covariate, length n
    vec beta   = as<vec>(largs["beta"]); // length 2
    int n      = time.size();
  
    vec grad(beta.size());
    vec risk;
  
    for (int i=0; i<n; ++i) {
      if (event(i) == 1) {
	risk = exp(beta(0)*x(span(i,n-1))+beta(1)*x(span(i,n-1))*log(time(i)));
	grad(0) += x(i) - sum(x(span(i,n-1)) % risk)/sum(risk);
	grad(1) += x(i)*log(time(i)) - sum(log(time(i))*x(span(i,n-1)) % risk)/sum(risk);
      }
    }
    return wrap(grad);
  }


  struct Tvc {
    int n;
    vec time, event, x, beta;
  };


  double test_cox_tvc3_negll(int n, double * beta, void * e) {
  
    Tvc * data = (Tvc *) e;

    double llike = 0.0;
  
    vec eta;
    // ivec index = linspace<ivec>(0, data->n - 1, data->n);
    // uvec indexi;

    for (int i=0; i < data->n; ++i) {
      if (data->event(i) == 1) {
	// indexi = find(index >= i);
	eta = beta[0]*data->x(span(i,data->n-1)) + 
	  beta[1]*log(data->time(i)) * data->x(span(i,data->n-1));
	llike += eta(0) - log(sum(exp(eta)));
      }
    }
    return -llike;
  }

  void test_cox_tvc3_negll_gr(int ncoef, double * beta, double * gr, void * e) {
  
    Tvc * data = (Tvc *) e;
  
    vec risk, subx;

    gr[0] = gr[1] = 0.0;
    
    for (int i=0; i<data->n; ++i) {
      if (data->event(i) == 1) {
	subx = data->x(span(i,data->n-1));
	risk = exp(beta[0]*subx + beta[1]*log(data->time(i)) * subx);
	gr[0] -= data->x(i) - sum(subx % risk)/sum(risk);
	gr[1] -= data->x(i)*log(data->time(i)) - sum(log(data->time(i)) * subx % risk)/sum(risk);
      }
    }
  }

  RcppExport SEXP test_cox_tvc3(SEXP args) {
  
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    vec x      = as<vec>(largs["x"]); // one covariate, length n
    NumericVector beta   = as<NumericVector>(largs["beta"]); // length 2
    int n      = time.size();

    Tvc data = {n, time, event, x, beta};

    BFGS bfgs;
    bfgs.optim(test_cox_tvc3_negll,test_cox_tvc3_negll_gr,beta, (void *) &data);

    return wrap(List::create(_("coef")=bfgs.coef,
			     _("negll")=bfgs.Fmin,
			     _("hessian")=bfgs.hessian));

  }

} // namespace
