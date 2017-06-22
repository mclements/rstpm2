#include <RcppArmadillo.h>
#include "c_optim.h"

namespace rstpm2 {

  using namespace Rcpp;
  using namespace arma;

  // assume right censored values in ascending order of time

  RcppExport SEXP test_cox_tvc2(SEXP args) {
  
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    mat X      = as<mat>(largs["X"]); // design matrix, n*c
    vec beta   = as<vec>(largs["beta"]); // length c+1
    int k      = as<int>(largs["k"]); // column to use for tvc
    int n      = time.size();
    int c      = X.n_cols;
    double llike = 0.0, lsum;
    vec eta;
    vec eta0 = X * beta(span(0,c-1));
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	eta = eta0(span(i,n-1)) + beta(c)*log(time(i))*X(span(i,n-1),k);
        lsum = log(sum(exp(eta)));
        for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
          llike += eta(j-i) - lsum;
        }
      }
    }
    return wrap(llike);
  }

  // Breslow estimator of the baseline cumulative hazard;
  // assumes right censored values in ascending order of time
  RcppExport SEXP test_cox_tvc2_H0(SEXP args) {
  
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    mat X      = as<mat>(largs["X"]); // design matrix, n*c
    vec beta   = as<vec>(largs["beta"]); // length c+1
    int k      = as<int>(largs["k"]); // column to use for tvc
    int n      = time.size();
    int c      = X.n_cols;
    double logHk, lsum;
    std::vector<double> H, etimes;
    vec eta;
    vec eta0 = X * beta(span(0,c-1));
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	eta = eta0(span(i,n-1)) + beta(c)*log(time(i))*X(span(i,n-1),k);
	logHk = 0.0;
        lsum = sum(exp(eta));
        for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
          logHk += 1 - lsum;
        }
	etimes.push_back(time(i));
	H.push_back(exp(logHk));
      }
    }
    return List::create(_("H")=wrap(H),_("times")=wrap(etimes));
  }


  // // Breslow estimator of the cumulative hazard for a given set of covariates;
  // // assumes right censored values in ascending order of time
  // RcppExport SEXP test_cox_tvc2_Hx0(SEXP args) {
  //   List largs = as<List>(args);
  //   vec time   = as<vec>(largs["time"]); // length n
  //   vec event  = as<vec>(largs["event"]); // length n
  //   mat X      = as<mat>(largs["X"]); // design matrix, n*c
  //   vec x0      = as<vec>(largs["x0"]); // design row for specific covariate pattern, length c
  //   vec beta   = as<vec>(largs["beta"]); // length c+1
  //   int k      = as<int>(largs["k"]); // column to use for tvc
  //   int n      = time.size();
  //   int c      = X.n_cols;
  //   double logHk, lsum;
  //   std::vector<double> H, etimes;
  //   vec eta;
  //   vec eta1 = X * beta(span(0,c-1));
  //   double eta0;
  //   for (int i=0; i<n; ++i) {
  //     if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
  // 	eta = eta1(span(i,n-1)) + beta(c)*log(time(i))*X(span(i,n-1),k);
  // 	eta0 = x0 * beta(span(0,c-1)) + beta(c)*log(time(i))*x0(k);
  // 	logHk = 0.0;
  //       lsum = sum(exp(eta));
  //       for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
  //         logHk += eta0 - lsum;
  //       }
  // 	etimes.push_back(time(i));
  // 	H.push_back(exp(logHk));
  //     }
  //   }
  //   return List::create(_("H")=wrap(H),_("times")=wrap(etimes));
  // }

  RcppExport SEXP test_cox_tvc2_grad(SEXP args) {
  
    List largs = as<List>(args);
    vec time   = as<vec>(largs["time"]); // length n
    vec event  = as<vec>(largs["event"]); // length n
    mat X      = as<mat>(largs["X"]); // one covariate, length n
    vec beta   = as<vec>(largs["beta"]); // length 2
    int k      = as<int>(largs["k"]); // column to use for tvc
    int n      = time.size();
    int c      = X.n_cols;
    vec grad(beta.size(),fill::zeros);
    vec lsum(beta.size(),fill::zeros);
    vec eta0 = X * beta(span(0,c-1));
    vec risk;
    mat Xrisk;
    for (int i=0; i<n; ++i) {
      if (event(i) == 1 && (i==0 || time(i) != time(i-1))) {
	risk = exp(eta0(span(i,n-1)) + beta(c)*log(time(i))*X(span(i,n-1),k));
        Xrisk = X(span(i,n-1), span::all);
        Xrisk.each_col() %= risk;
        lsum(span(0,c-1)) = sum(Xrisk, 0).t() / sum(risk);
        lsum(c) = sum(log(time(i))*Xrisk(span::all,k))/sum(risk);
        for (int j=i; j<n && time(j)==time(i) && event(j)==1; ++j) {
	  grad(span(0,c-1)) += X(j,span::all).t() - lsum(span(0,c-1));
	  grad(c) += X(j,k)*log(time(i)) - lsum(c);
        }
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

    for (int i=0; i < data->n; ++i) {
      if (data->event(i) == 1) {
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
    NumericVector beta   = clone(as<NumericVector>(largs["beta"])); // length 2
    int n      = time.size();

    Tvc data = {n, time, event, x, beta};

    BFGS bfgs;
    bfgs.optim(test_cox_tvc3_negll,test_cox_tvc3_negll_gr,beta, (void *) &data);

    return wrap(List::create(_("coef")=bfgs.coef,
			     _("negll")=bfgs.Fmin,
			     _("hessian")=bfgs.hessian));

  }

} // namespace
