#include <RcppArmadillo.h>
#include <vector>
#include <map>
#include "c_optim.h"

namespace rstpm2 {

  // import namespaces
  using namespace Rcpp;
  using namespace arma;
  // typedefs
  typedef bool constraintfn(int, double *, void *);
  // enums
  enum IntervalCensoring {RightCensored, ExactTime, LeftCensored, IntervalCensored};
  // structs
  struct li_constraint { vec li; double constraint; };
  struct gradli_constraint { mat gradli; mat constraint; };
  li_constraint operator+(li_constraint const &left, li_constraint const &right) {
    li_constraint out = {left.li+right.li,left.constraint+right.constraint};
    return out;
  }
  gradli_constraint operator+(gradli_constraint const &left, gradli_constraint const &right) {
    gradli_constraint out = {left.gradli+right.gradli,left.constraint+right.constraint};
    return out;
  }
  // Hadamard element-wise multiplication for the _columns_ of a matrix with a vector
  mat rmult(mat m, vec v) {
    mat out(m);
    out.each_col() %= v;
    return out;
  }
  mat lmult(vec v, mat m) {
    return rmult(m,v);
  }
  // print utilities
  void Rprint(NumericMatrix const & m) {
    for (int i=0; i<m.nrow(); ++i) {
      for (int j=0; j<m.ncol(); ++j) 
	Rprintf("%f ", m(i,j));
      Rprintf("\n");
    }
  }
  void Rprint(NumericVector const & v) {
    for (int i=0; i<v.size(); ++i) 
      Rprintf("%f ", v(i));
    Rprintf("\n");
  }
  void Rprint(vec const & v) {
    for (size_t i=0; i<v.size(); ++i) 
      Rprintf("%f ", v(i));
    Rprintf("\n");
  }
  void Rprint(rowvec const & v) {
    for (size_t i=0; i<v.size(); ++i) 
      Rprintf("%f ", v(i));
    Rprintf("\n");
  }
  void Rprint(mat const & m) {
    for (size_t i=0; i<m.n_rows; ++i) {
      for (size_t j=0; j<m.n_cols; ++j) 
	Rprintf("%f ", m(i,j));
      Rprintf("\n");
    }
  }
  // vectorised functions
  vec pnorm01(vec const & x) {
    vec out(x.size());
    for (size_t i=0; i<x.size(); ++i)
      out(i) = R::pnorm(x(i),0.0,1.0,1,0);
    return out;
  }
  vec qnorm01(vec const & x) {
    vec out(x.size());
    for (size_t i=0; i<x.size(); ++i)
      out(i) = R::qnorm(x(i),0.0,1.0,1,0);
    return out;
  }
  vec dnorm01(vec const & x) {
    vec out(x.size());
    for (size_t i=0; i<x.size(); ++i)
      out(i) = R::dnorm(x(i),0.0,1.0,0);
    return out;
  }
  // we could use templates for the following...
  vec logit(vec const & p) {
    return log(p/(1-p));
  }
  vec expit(vec const & x) {
    return 1/(1+exp(-x));
  }
  double logit(double p) {
    return log(p/(1-p));
  }
  double expit(double x) {
    return 1/(1+exp(-x));
  }
  // utility for fast ragged sums
  RcppExport SEXP tapplySum(SEXP s_y, SEXP s_group) {
    NumericVector y(s_y);
    NumericVector group(s_group);
    NumericVector::iterator iy, igroup;
    std::map<double,double> out;
    for (iy = y.begin(), igroup = group.begin(); iy != y.end(); ++iy, ++igroup) { 
	out[*igroup] += *iy;
    }
    return wrap(out);
  }

  // Various link objects
  class LogLogLink {
  public:
    vec link(vec S) { return log(-log(S)); }
    vec ilink(vec x) { return exp(-exp(x)); }
    vec h(vec eta, vec etaD) { return etaD % exp(eta); }
    vec H(vec eta) { return exp(eta); }
    mat gradH(vec eta, mat X) { return rmult(X,exp(eta)); }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(XD, exp(eta)) + rmult(X, etaD % exp(eta));
    }
    cube hessianH(vec beta, mat X) { 
      cube c(beta.size(), beta.size(), X.n_rows);
      for (size_t i=0; i<X.n_rows; ++i) 
	c.slice(i) = X.row(i).t()*X.row(i)*exp(dot(X.row(i),beta)); 
      return c;
    }
    cube hessianh(vec beta, mat X, mat XD) {
      cube c(beta.size(), beta.size(), X.n_rows);
      for (size_t i=0; i<X.n_rows; ++i) {
	rowvec Xi = X.row(i);
	rowvec XDi = XD.row(i);
	c.slice(i) = (XDi.t()*Xi + Xi.t()*XDi + dot(XDi,beta)*Xi.t()*Xi)*exp(dot(Xi,beta));
      }
      return c;
    }
  };
  // Useful relationship: d/dx expit(x)=expit(x)*expit(-x) 
  class LogitLink {
  public:
    vec link(vec S) { return -logit(S); }
    vec ilink(vec x) { return expit(-x); }
    // vec h(vec eta, vec etaD) { return etaD % exp(eta) % expit(-eta); }
    vec h(vec eta, vec etaD) { return etaD % expit(eta); }
    vec H(vec eta) { return -log(expit(-eta)); }
    mat gradH(vec eta, mat X) { 
      return rmult(X,expit(eta));
    }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      // return rmult(X, etaD % exp(eta) % expit(-eta)) - 
      // 	rmult(X, exp(2*eta) % etaD % expit(-eta) % expit(-eta)) +
      // 	rmult(XD, exp(eta) % expit(-eta));
      return rmult(XD, expit(eta)) + 
	rmult(X, expit(eta) % expit(-eta) % etaD);
    }
    cube hessianH(vec beta, mat X) { 
      cube c(beta.size(), beta.size(), X.n_rows);
      for (size_t i=0; i<X.n_rows; ++i) {
	colvec Xi = X.row(i);
	double Xibeta = dot(Xi,beta);
	c.slice(i) = Xi.t()*Xi*expit(Xibeta)*expit(-Xibeta); 
      }
      return c;
    }
    cube hessianh(vec beta, mat X, mat XD) { 
      cube c(beta.size(), beta.size(), X.n_rows);
      for (size_t i=0; i<X.n_rows; ++i) {
	colvec Xi = X.row(i);
	colvec XDi = XD.row(i);
	double Xibeta = dot(Xi,beta);
	c.slice(i) = XDi.t()*Xi*expit(Xibeta) + Xi.t()*XDi*expit(-Xibeta)*expit(Xibeta); 
      }
      return c;
    }
  };
  class ProbitLink {
  public:
    vec link(vec S) { return -qnorm01(S); }
    vec ilink(vec eta) { return pnorm01(-eta); }
    vec H(vec eta) { return -log(pnorm01(-eta)); }
    vec h(vec eta, vec etaD) { return etaD % dnorm01(-eta) / pnorm01(-eta); }
    mat gradH(vec eta, mat X) { 
      return rmult(X, dnorm01(-eta) / pnorm01(-eta));
    }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(X, -eta % dnorm01(eta) % etaD / pnorm01(-eta)) +
	rmult(X, dnorm01(eta) % dnorm01(eta) / pnorm01(-eta) / pnorm01(-eta) % etaD) +
	rmult(XD,dnorm01(eta) / pnorm01(-eta));
    }
  };
  class LogLink {
  public:
    vec link(vec S) { return -log(S); }
    vec ilink(vec x) { return exp(-x); }
    vec h(vec eta, vec etaD) { return etaD; }
    vec H(vec eta) { return eta; }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { return XD; }
    mat gradH(vec eta, mat X) { return X;  }
  };

  // wrappers to call class methods from C
  template<class T>
  double optimfunction(int n, double * beta, void * void_obj) {
    T * obj = static_cast<T*>(void_obj);
    vec coef(beta,n);
    double value = obj->objective(coef % obj->parscale);
    if (obj->bfgs.trace>1) {
      Rprintf("beta="); Rprint(coef);
      Rprintf("objective=%g\n",value);
    };
    return value;
  }
  template<class T>
  void optimgradient(int n, double * beta, double * grad, void * void_obj) {
    T * obj = static_cast<T*>(void_obj);
    vec coef(beta,n);
    if (obj->bfgs.trace>1) {
      Rprintf("beta="); Rprint(coef);
    }
    if (obj->bfgs.trace>2) {
      Rprintf("parscale="); Rprint(obj->parscale);
    }
    vec gr = obj->gradient(coef % obj->parscale);
    if (obj->bfgs.trace>1) {
      Rprintf("gradient="); Rprint(gr);
    }
    if (obj->bfgs.trace>2) {
      Rprintf("fdgradient="); Rprint(obj->fdgradient(coef % obj->parscale));
    }
    for (int i=0; i<n; ++i) grad[i] = gr[i];
  }
  template<class T>
  void optimfunction_nlm(int n, double * beta, double * f, void * void_obj) {
    T * obj = static_cast<T*>(void_obj);
    vec coef(beta,n);
    *f = obj->objective(coef % obj->parscale);
  }
  class BFGS2 : public BFGS {
  public:
    NumericMatrix calc_hessian(optimgr gr, void * ex) {
      if (parscale.size()==0) REprintf("parscale is not defined for BFGS2::calc_hessian.\n");
      int n = coef.size();
      NumericVector df1(n);
      NumericVector df2(n);
      NumericMatrix hess(n,n);
      double tmp;
      for(int i=0; i<n; ++i) {
	tmp = coef[i];
	coef[i] = tmp + epshess/parscale[i];
	gr(n, &coef[0], &df1[0], ex);
	coef[i] = tmp - epshess/parscale[i];
	gr(n, &coef[0], &df2[0], ex);
	for (int j=0; j<n; ++j)
	  hess(i,j) = (df1[j] - df2[j]) / (2*epshess);
	coef[i] = tmp;
      }
      // now symmetrize
      for(int i=0; i<n; ++i) 
	for(int j=i; j<n; ++j) 
	  if (i != j)
	    hess(i,j) = hess(j,i) = (hess(i,j) + hess(j,i)) / 2.0;
      return hess; // wrap()?
    }
    vec parscale;
  };
  class NelderMead2 : public NelderMead {
  public:
    NumericMatrix calc_hessian(optimfn fn, void * ex) {
      if (parscale.size()==0) REprintf("parscale is not defined for NelderMead2::calc_hessian.");
      int n = coef.size();
      NumericMatrix hess(n,n);
      double tmpi,tmpj,f1,f0,fm1,hi,hj,fij,fimj,fmij,fmimj;
      f0 = fn(n,&coef[0],ex);
      for(int i=0; i<n; ++i) {
	tmpi = coef[i];
	hi = epshess*(1.0+std::abs(tmpi))/parscale[i];
	coef[i] = tmpi + hi;
	f1=fn(n, &coef[0], ex);
	coef[i] = tmpi - hi;
	fm1=fn(n, &coef[0], ex);
	// hess(i,i) = (-f2 +16.0*f1 - 30.0*f0 + 16.0*fm1 -fm2)/(12.0*hi*hi);
	hess(i,i) = (f1 - 2.0*f0 + fm1)/(hi*hi*parscale[i]*parscale[i]);
	coef[i] = tmpi;
	for (int j=i; j<n; ++j) {
	  if (i != j) {
	    tmpj = coef[j];
	    hj = epshess*(1.0+std::abs(tmpj))/parscale[j];
	    coef[i] = tmpi + hi;
	    coef[j] = tmpj + hj;
	    fij=fn(n, &coef[0], ex);
	    coef[i] = tmpi + hi;
	    coef[j] = tmpj - hj;
	    fimj=fn(n, &coef[0], ex);
	    coef[i] = tmpi - hi;
	    coef[j] = tmpj + hj;
	    fmij=fn(n, &coef[0], ex);
	    coef[i] = tmpi - hi;
	    coef[j] = tmpj - hj;
	    fmimj=fn(n, &coef[0], ex);
	    hess(j,i) = hess(i,j) = (fij-fimj-fmij+fmimj)/(4.0*hi*hj*parscale[i]*parscale[j]);
	    coef[i] = tmpi;
	    coef[j] = tmpj;
	  }
	}
      }
      return hess;
    }
    vec parscale;
  };

  class Nlm2 : public Nlm {
  public:
    NumericMatrix calc_hessian(fcn_p fn, void * ex) {
      if (parscale.size()==0) REprintf("parscale is not defined for Nlm2::calc_hessian.");
      int n = coef.size();
      NumericMatrix hess(n,n);
      double tmpi,tmpj,f1,f0,fm1,hi,hj,fij,fimj,fmij,fmimj;
      fn(n,&coef[0],&f0,ex);
      for(int i=0; i<n; ++i) {
  tmpi = coef[i];
	hi = epshess*(1.0+std::abs(tmpi))/parscale[i];
	coef[i] = tmpi + hi;
	fn(n, &coef[0], &f1, ex);
	coef[i] = tmpi - hi;
	fn(n, &coef[0], &fm1, ex);
	// hess(i,i) = (-f2 +16.0*f1 - 30.0*f0 + 16.0*fm1 -fm2)/(12.0*hi*hi);
	hess(i,i) = (f1 - 2.0*f0 + fm1)/(hi*hi*parscale[i]*parscale[i]);
	coef[i] = tmpi;
	for (int j=i; j<n; ++j) {
	  if (i != j) {
	    tmpj = coef[j];
	    hj = epshess*(1.0+std::abs(tmpj))/parscale[j];
	    coef[i] = tmpi + hi;
	    coef[j] = tmpj + hj;
	    fn(n, &coef[0], &fij, ex);
	    coef[i] = tmpi + hi;
	    coef[j] = tmpj - hj;
	    fn(n, &coef[0], &fimj, ex);
	    coef[i] = tmpi - hi;
	    coef[j] = tmpj + hj;
	    fn(n, &coef[0], &fmij, ex);
	    coef[i] = tmpi - hi;
	    coef[j] = tmpj - hj;
	    fn(n, &coef[0], &fmimj, ex);
	    hess(j,i) = hess(i,j) = (fij-fimj-fmij+fmimj)/(4.0*hi*hj*parscale[i]*parscale[j]);
	    coef[i] = tmpi;
	    coef[j] = tmpj;
	  }
	}
      }
      return hess;
    }
    vec parscale;
  };


  // Main class for Stpm2 (fully parametric)
  template<class Link = LogLogLink>
  class Stpm2 : public Link {
  public:
    typedef Stpm2<Link> This;
    using Link::h;
    using Link::H;
    using Link::gradh;
    using Link::gradH;
    Stpm2(SEXP sexp) : bfgs() {
      List list = as<List>(sexp);
      bfgs.coef = init = as<NumericVector>(list["init"]);
      X = as<mat>(list["X"]); 
      XD = as<mat>(list["XD"]); 
      bhazard = as<vec>(list["bhazard"]);
      wt = as<vec>(list["wt"]);
      event = as<vec>(list["event"]);
      time = as<vec>(list["time"]);
      delayed = as<bool>(list["delayed"]);
      interval = as<bool>(list["interval"]);
      n = init.size(); // number of parameters
      N = X.n_rows; // number of observations
      X1 = as<mat>(list["X1"]);
      X0 = as<mat>(list["X1"]);
      wt0.set_size(N); wt0.fill(0.0);
      if (delayed) {
	X0 = as<mat>(list["X0"]);
	wt0 = as<vec>(list["wt0"]);
      }
      if (interval) {
	X1 = as<mat>(list["X1"]);
      }
      ind0 = as<uvec>(list["ind0"]); // length N boolean
      map0 = as<uvec>(list["map0"]); // length N map from individuals to row of X0
      which0 = as<uvec>(list["which0"]); // length N0 indicator for X0
      kappa = as<double>(list["kappa"]);
      optimiser = as<std::string>(list["optimiser"]);
      bfgs.trace = as<int>(list["trace"]);
      reltol = bfgs.reltol = as<double>(list["reltol"]);
      bfgs.hessianp = false;
      bfgs.parscale = parscale = as<vec>(list["parscale"]);
    }
    // log-likelihood components and constraints
    // Note: the li and gradli components depend on other variables (bhazard, wt, wt0, wt, event);
    //       calculations should *not* be done on subsets
    li_constraint li_right_censored(vec eta, vec etaD) {
      vec h = Link::h(eta,etaD) + bhazard;
      vec H = Link::H(eta);
      vec eps = h*0.0 + 1.0e-16; 
      double constraint = kappa/2.0 * (sum(h % h % (h<0)) +
				       sum(H % H % (H<0))); 
      h = max(h,eps);
      H = max(H,eps);
      vec li = wt % event % log(h) - wt % H;
      li_constraint out = {li, constraint};
      return out;
    }
    li_constraint li_left_truncated(vec eta0) {
      vec H0 = Link::H(eta0);
      double constraint = kappa/2.0 * sum(H0 % H0 % (H0<0));
      vec eps = H0*0.0 + 1e-16;
      vec li = wt0 % max(H0,eps);
      li_constraint out = {li, constraint};
      return out;
    }
    li_constraint li_interval(vec eta, vec etaD, vec eta1) {
      vec H = Link::H(eta);
      vec H1 = Link::H(eta1);
      vec h = Link::h(eta, etaD) + bhazard;
      double constraint = kappa/2.0 * (sum(H % H % (H<0))+
				       sum(h % h % (h<0))+
				       sum(H1 % H1 % (H1<0)));
      vec eps = H*0.0 + 1e-16;
      H = max(H,eps);
      H1 = max(H1,eps);
      h = max(h,eps);
      vec li(N,fill::zeros);
      uvec index;
      index = find(ttype == 0); // right censored
      if (any(index)) li(index) -= wt(index) % H(index);
      index = find(ttype == 1); // exact
      if (any(index)) li(index) += wt(index) % (log(h(index)) - H(index));
      index = find(ttype == 2); // left censored
      if (any(index)) li(index) += wt(index) % log(1-exp(-H(index)));
      index = find(ttype == 3); // interval censored
      if (any(index)) li(index) += wt(index) % log(exp(-H(index)) - exp(-H1(index)));
      li_constraint out = {li, constraint};
      return out;
    }
    li_constraint li(vec eta, vec etaD, vec eta0, vec eta1) {
      if (interval) {
	return li_interval(eta, etaD, eta1);
      }
      else {
	li_constraint s = li_right_censored(eta, etaD);
	if (delayed) {
	  li_constraint s0 = li_left_truncated(eta0);
	  s.constraint += s0.constraint;
	  s.li(which0) += s0.li;
	}
	return s;
      }
    }
    // negative log-likelihood
    double objective(vec beta) {
      li_constraint s = li(X * beta, XD * beta, X0 * beta, X1 * beta);
      return -sum(s.li) + s.constraint;
    }
    // finite-differencing of the gradient for the objective
    vec fdgradient(vec beta, double eps = 1.0e-8) {
      int n = beta.size();
      vec grad(n);
      for (int i=0; i<n; ++i) {
	vec lower(beta);
	vec upper(beta);
	upper(i) += eps;
	lower(i) -= eps;
	grad(i) = (objective(upper) - objective(lower)) / 2.0 / eps;
      }
      return grad;
    }
    // log-likelihood gradient components
    // Note: gradH and gradh are currently calculated at eta and etaD and may not satisfy the constraints, 
    //       while H and h may be transformed for the constraints 
    gradli_constraint gradli_right_censored(vec eta, vec etaD, mat X, mat XD) {
      vec h = Link::h(eta,etaD) + bhazard;
      vec H = Link::H(eta);
      vec eps = h*0.0 + 1.0e-16; // hack
      mat gradH = Link::gradH(eta,X);
      mat gradh = Link::gradh(eta,etaD,X,XD);
      mat Xconstraint = kappa * (rmult(gradh,h%(h<eps)) +
				 rmult(gradH,H % (H<eps)));
      H = max(H,eps);
      h = max(h,eps);
      mat Xgrad = -rmult(gradH, wt) + rmult(gradh, event / h % wt);
      gradli_constraint out = {Xgrad, Xconstraint};
      return out;
    }
    gradli_constraint gradli_left_truncated(vec eta0, mat X0) {
      mat gradH0 = Link::gradH(eta0, X0); 
      vec H0 = Link::H(eta0); 
      vec eps = H0*0.0 + 1.0e-16; // hack
      mat Xconstraint = kappa * rmult(gradH0, H0 % (H0<eps));
      H0 = max(H0,eps);
      mat Xgrad = rmult(gradH0, wt0);
      gradli_constraint out = {Xgrad, Xconstraint};
      return out;
    }
    gradli_constraint gradli_interval_censored(vec eta, vec etaD, vec eta1, 
						       mat X, mat XD, mat X1) {
      vec H = Link::H(eta);
      vec h = Link::h(eta,etaD);
      vec H1 = Link::H(eta1);
      mat gradH = Link::gradH(eta,X);
      mat gradH1 = Link::gradH(eta1,X1);
      mat gradh = Link::gradh(eta, etaD, X, XD);
      vec eps = H*0.0 + 1e-16;
      mat Xconstraint = kappa * (rmult(gradH, H % (H<eps))+
				 rmult(gradh, h % (h<eps))+
				 rmult(gradH1, H1 % (H1<eps)));
      H = max(H,eps);
      H1 = max(H1,eps);
      h = max(h,eps);
      mat li(N,n,fill::zeros);
      uvec index;
      index = find(ttype == 0); // right censored
      if (any(index)) li(index) -= rmult(gradH.rows(index),wt(index));
      index = find(ttype == 1); // exact
      if (any(index)) li(index) += rmult(gradh.rows(index),wt(index) / h(index)) - rmult(gradH.rows(index),wt(index));
      index = find(ttype == 2); // left censored
      if (any(index)) li(index) += rmult(-gradH.rows(index),-exp(-H(index)) / (1-exp(-H(index))) % wt(index));
      index = find(ttype == 3); // interval censored
      if (any(index)) {
	vec V = wt(index) / (exp(-H(index)) - exp(-H1(index)));
	li(index) += rmult(gradH1.rows(index),V % exp(-H1(index))) - 
	  rmult(gradH.rows(index),V % exp(-H(index)));
      }
      gradli_constraint out = {li, Xconstraint};
      return out;
    }
    gradli_constraint gradli(vec eta, vec etaD, vec eta0, vec eta1,
				     mat X, mat XD, mat X0, mat X1) {
      if (interval) return gradli_interval_censored(eta, etaD, eta1, X, XD, X1);
      else {
	gradli_constraint s = gradli_right_censored(eta, etaD, X, XD);
	if (delayed) {
	  gradli_constraint s0 = gradli_left_truncated(eta0, X0);
	  s.constraint.rows(which0) += s0.constraint;
	  s.gradli.rows(which0) += s0.gradli;
	}
	return s;
      }
    }
    // gradient of the negative log-likelihood
    vec gradient(vec beta) {
      gradli_constraint gc = gradli(X * beta, XD * beta, X0 * beta, X1 * beta,
				    X, XD, X0, X1);
      rowvec dconstraint = sum(gc.constraint,0);
      rowvec vgr = sum(gc.gradli,0);
      vec gr(n);
      for (size_t i = 0; i<beta.size(); ++i) {
	gr[i] = vgr[i] - dconstraint[i];
      }
      return -gr;
    }
    bool feasible(vec beta) {
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec h = Link::h(eta, etaD) + bhazard;
      vec H = Link::H(eta);
      bool condition = all((h>0) % (H>0));
      if (delayed) {
	vec eta0 = X0 * beta;
	// vec etaD0 = XD0 * beta;
	vec H0 = Link::H(eta0);
	condition = condition && all(H0>0);
      }
      return condition;
    }
    void pre_process() {
      for (int i = 0; i<n; ++i) init[i] /= parscale[i];
    }
    void post_process() {
      for (int i = 0; i<n; ++i) {
	bfgs.coef[i] *= parscale[i];
	init[i] *= parscale[i];
      }
    }
    void optimWithConstraintBFGS(NumericVector init) {
      bool satisfied;
      do {
	bfgs.optim(&optimfunction<This>, &optimgradient<This>, init, (void *) this);
	vec vcoef(&bfgs.coef[0],n);
	satisfied = feasible(vcoef % parscale);
	if (!satisfied) kappa *= 2.0;   
      } while ((!satisfied) && kappa < 1.0e3);
    }
    void optimWithConstraintNM(NumericVector init) {
      bool satisfied;
      NelderMead2 nm;
      nm.hessianp = false;
      nm.parscale = parscale;
      do {
	nm.optim(&optimfunction<This>, init, (void *) this);
	vec vcoef(&nm.coef[0],n);
	satisfied = feasible(vcoef % parscale);
	if (!satisfied) kappa *= 2.0;   
      } while ((!satisfied) && kappa < 1.0e3);
      nm.hessian = nm.calc_hessian(&optimfunction<This>, (void *) this);
      bfgs.coef = nm.coef;
      bfgs.hessian = nm.hessian;
    }
    void optimWithConstraint(NumericVector init) {
      if (this->optimiser == "NelderMead")
	optimWithConstraintNM(init);
      else
	optimWithConstraintBFGS(init);
    }
    // find the left truncated values for a given index
    uvec find0() { return map0(find(ind0)); }
    uvec find0(vec index) { return (map0(index))(find(ind0(index))); }
    NumericVector init;
    mat X, XD, X0, X1; 
    vec bhazard,wt,wt0,event,time,parscale,ttype;
    double kappa, reltol;
    bool delayed, interval;
    int n, N;
    BFGS2 bfgs;
    uvec ind0, map0, which0;
    std::string optimiser;
  };

  // RcppExport SEXP optim_stpm2_nlm(SEXP args) {
  //   stpm2 data(args);
  //   data.init = ew_div(data.init,data.parscale);
  //   int n = data.init.size();
  //   Nlm nlm;
  //   nlm.gradtl = nlm.steptl = data.reltol;
  //   //nlm.method=2; nlm.stepmx=0.0;
  //   bool satisfied;
  //   do {
  //     nlm.optim(& fminfn_nlm<stpm2>, & grfn<stpm2>, data.init, (void *) &data);
  //     satisfied = fminfn_constraint<stpm2>(n,&nlm.coef[0],(void *) &data);
  //     if (!satisfied) data.kappa *= 2.0;
  //   } while (!satisfied && data.kappa<1.0e5);
  //   nlm.coef = ew_mult(nlm.coef, data.parscale);
  //   return List::create(_("itrmcd")=wrap(nlm.itrmcd),
  // 			_("niter")=wrap(nlm.itncnt),
  // 			_("coef")=wrap(nlm.coef),
  // 			_("hessian")=wrap(nlm.hessian));
  // }

  // penalised smoothers 
  // _Not_ sub-classes of Pstpm2, hence the longer function signatures
  class SmoothLogH {
  public:
    struct Smoother {
      int first_para, last_para;
      mat S;
    };
    SmoothLogH(SEXP sexp) {
      List list = as<List>(sexp);
      List lsmooth = as<List>(list["smooth"]);
      for(int i=0; i<lsmooth.size(); ++i) {
	List lsmoothi = as<List>(lsmooth[i]);
	List lsmoothiS = as<List>(lsmoothi["S"]);
	Smoother smoothi =  {
	  as<int>(lsmoothi["first.para"]) - 1, 
	  as<int>(lsmoothi["last.para"]) - 1, 
	  as<mat>(lsmoothiS[0])
	};
	smooth.push_back(smoothi);
      }
    }
    double penalty(vec vbeta, vec sp) {
      double value = 0.0;
      for (size_t i=0; i < smooth.size(); ++i) {
	Smoother smoothi = smooth[i];
	value += sp[i]/2 * 
	  dot(vbeta.subvec(smoothi.first_para,smoothi.last_para),
	      smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para));
      }
      return value;  
    }
    vec penalty_gradient(vec vbeta, vec sp) {
      int n = vbeta.size();
      rowvec vgr(n, fill::zeros);
      for (size_t i=0; i < smooth.size(); ++i) {
	Smoother smoothi = smooth[i];
	vgr.subvec(smoothi.first_para,smoothi.last_para) += 
	  sp[i] * (smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para)).t();
      }
      vec gr(n);
      for (int i=0; i<n; ++i) gr[i] = vgr[i];
      return gr;
    }
    vec traces(mat X) {
      vec values(smooth.size(), fill::zeros);
      for (size_t i=0; i < smooth.size(); ++i) {
	Smoother smoothi = smooth[i];
	for (size_t j=smoothi.first_para; j <= static_cast<size_t>(smoothi.last_para); ++j)
	  values[i] += X(j,j);
      }
      return values;
    }
    std::vector<Smoother> smooth;
  };
  // // TODO: include first_para and last_para class data for the traces() method
  // class SmoothHaz {
  // public:
  //   struct Smoother {
  //     mat X0, X1, X2, X3;
  //     vec w;
  //     double lambda;
  //   };
  //   SmoothHaz(SEXP sexp) {
  //     List list = as<List>(sexp);
  //     List lsmooth = as<List>(list["smooth"]);
  //     for(int i=0; i<lsmooth.size(); ++i) {
  // 	List lsmoothi = as<List>(lsmooth[i]);
  // 	Smoother smoothi =  {
  // 	  as<mat>(lsmoothi["X0"]), 
  // 	  as<mat>(lsmoothi["X1"]), 
  // 	  as<mat>(lsmoothi["X2"]), 
  // 	  as<mat>(lsmoothi["X3"]), 
  // 	  as<vec>(lsmoothi["w"]), 
  // 	  as<double>(lsmoothi["lambda"])
  // 	};
  // 	smooth.push_back(smoothi);
  //     }
  //   }
  //   double penalty(vec vbeta, vec sp) {
  //     double value = 0.0;
  //     for (size_t i=0; i < smooth.size(); ++i) {
  // 	Smoother obj = smooth[i];
  // 	vec s0 = obj.X0 * vbeta;
  // 	vec s1 = obj.X1 * vbeta;
  // 	vec s2 = obj.X2 * vbeta;
  // 	vec s3 = obj.X3 * vbeta;
  // 	vec h2 = exp(s0) % (s3 + 3*s1%s2 + s1%s1%s1);
  // 	value += sp[i]/2 * obj.lambda * sum(obj.w % h2 % h2);
  //     }
  //     return value;
  //   }
  // vec penalty_gradient(vec vbeta, vec sp) {
  //   int n = vbeta.size();
  //   rowvec vgr(n, fill::zeros);
  //   for (size_t i=0; i < smooth.size(); ++i) {
  //     Smoother obj = smooth[i];
  //     vec s0 = obj.X0 * vbeta;
  //     vec s1 = obj.X1 * vbeta;
  //     vec s2 = obj.X2 * vbeta;
  //     vec s3 = obj.X3 * vbeta;
  //     vec h2 = exp(s0) % (s3 + 3*s1%s2 + s1%s1%s1);
  //     mat dh2sq_dbeta = 
  // 	2*lmult(h2 % exp(s0),
  // 		obj.X3+3*(rmult(obj.X1,s2)+rmult(obj.X2,s1))+3*lmult(s1%s1,obj.X1))+
  // 	2*lmult(h2 % h2,obj.X0);
  //     vgr += sp[i]*obj.lambda*sum(lmult(obj.w, dh2sq_dbeta),0);
  //   }
  //   vec gr(n);
  //   for (int i = 0; i<n; ++i) gr[i] = vgr[i];
  //   return gr;
  // }
  //   std::vector<Smoother> smooth;
  // };

  template<class T>
  double pstpm2_multivariate_step(int n, double * logsp_ptr, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    vec logsp(logsp_ptr,n);
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
    return model->multivariate_step(logsp);
  }    
  template<class T>
  double pstpm2_first_step(double logsp, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
    return model->first_step(logsp);
  }    
  template<class T>
  void pstpm2_multivariate_stepNlm(int n, double * logsp_ptr, double * objective, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    vec logsp(logsp_ptr,n);
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
    *objective = model->multivariate_step(logsp);
  }    

  // Penalised link-based models
  template<class Stpm2Type = Stpm2<LogLogLink>, class Smooth = SmoothLogH>
  class Pstpm2 : public Stpm2Type, public Smooth {
  public:
    typedef Pstpm2<Stpm2Type,Smooth> This;
    Pstpm2(SEXP sexp) : Stpm2Type(sexp), Smooth(sexp) {
      List list = as<List>(sexp);
      sp = as<vec>(list["sp"]);
      reltol_search = as<double>(list["reltol_search"]);
      reltol_outer = as<double>(list["reltol_outer"]);
      alpha = as<double>(list["alpha"]);
      criterion = as<int>(list["criterion"]);
      outer_optim = as<int>(list["outer_optim"]);
    }
    double objective(vec beta) {
      return Stpm2Type::objective(beta) + Smooth::penalty(beta,sp);
    }
    double objective0(vec beta) {
      return Stpm2Type::objective(beta);
    }
    vec gradient(vec beta) {
      return Stpm2Type::gradient(beta) + Smooth::penalty_gradient(beta,sp);
    }
    vec gradient0(vec beta) {
      return Stpm2Type::gradient(beta);
    }
    // is the following strictly needed - or even correct?
    void optimWithConstraint(NumericVector init) {
      bool satisfied;
      if (this->bfgs.trace > 0) 
	Rprintf("Starting optimization\n");
      do {
	this->bfgs.optim(&optimfunction<This>, &optimgradient<This>, init, (void *) this);
	vec vcoef(&this->bfgs.coef[0],this->n);
	satisfied = Stpm2Type::feasible(vcoef % this->parscale);
	if (!satisfied) this->kappa *= 2.0;   
      } while ((!satisfied) && this->kappa < 1.0e3);
    }
    double first_step(double logsp) {
      sp[0] = exp(logsp);
      this->pre_process();
      this->optimWithConstraint(this->init);
      this->bfgs.hessian = this->bfgs.calc_hessian(&optimgradient<This>, (void *) this);
      NumericMatrix hessian0 = this->bfgs.calc_hessian(&optimgradient<Stpm2Type>, (void *) this);
      if (this->bfgs.trace > 1)  {
	Rprintf("Debug on trace calculation. Coef:\n");
	Rprint(this->bfgs.coef);
	if (this->bfgs.trace > 1) {
	  Rprintf("Hessian0:\n");
	  Rprint(hessian0);
	  Rprintf("Hessian:\n");
	  Rprint(this->bfgs.hessian);
	}
      }
      double edf = arma::trace(solve(as<mat>(this->bfgs.hessian),as<mat>(hessian0)));
      double negll = this->bfgs.calc_objective(&optimfunction<Stpm2Type>, (void *) this);
      double gcv =  negll + alpha*edf;
      double bic =  negll + alpha*edf*log(sum(this->event));
      this->init = this->bfgs.coef;
      if (this->bfgs.trace > 0)
	Rprintf("sp=%f\tedf=%f\tnegll=%f\tgcv=%f\tbic=%f\talpha=%f\n",sp[0],edf,negll,gcv,bic,alpha);
      this->post_process();
      return criterion==1 ? gcv : bic;
    }
    double multivariate_step(vec logsp) {
      this->pre_process();
      double lsp = -9.0, usp = 9.0;
      for (size_t i=0; i < sp.size(); ++i)
	sp[i] = exp(bound(logsp[i],lsp,usp));
      if (this->bfgs.trace > 0) {
	for (size_t i = 0; i < sp.size(); ++i)
	  Rprintf("sp[%i]=%f\t",i,sp[i]);
      }
      this->optimWithConstraint(this->init);
      this->bfgs.hessian = this->bfgs.calc_hessian(&optimgradient<This>, (void *) this);
      NumericMatrix hessian0 = this->bfgs.calc_hessian(&optimgradient<Stpm2Type>, (void *) this);
      double edf = arma::trace(solve(as<mat>(this->bfgs.hessian),as<mat>(hessian0)));
      double negll = this->bfgs.calc_objective(&optimfunction<Stpm2Type>, (void *) this);
      double gcv =  negll + alpha*edf;
      double bic =  2.0*negll + alpha*edf*log(sum(this->event));
      this->init = this->bfgs.coef;
      // simple boundary constraints
      double constraint = 0.0;
      // double kappa = 1.0;
      for (size_t i=0; i < sp.size(); ++i) {
	if (logsp[i] < lsp)
	  constraint +=  pow(logsp[i]-lsp,2.0);
	if (logsp[i] > usp)
	  constraint +=  pow(logsp[i]-usp,2.0);
      }
      double objective =  criterion==1 ? 
	gcv + this->kappa/2.0*constraint : 
	bic + this->kappa/2.0*constraint;
      if (this->bfgs.trace > 0) {
	Rprintf("edf=%f\tnegll=%f\tgcv=%f\tbic=%f\tobjective=%f\n",edf,negll,gcv,bic,objective);
      }
      this->post_process();
      return objective;
    }
    SEXP optim_fixed() { 
      this->pre_process();
      this->optimWithConstraint(this->init);
      this->bfgs.hessian = this->bfgs.calc_hessian(&optimgradient<This>, (void *) this);
      NumericMatrix hessian0 = this->bfgs.calc_hessian(&optimgradient<Stpm2Type>, (void *) this);
      mat Proj = solve(as<mat>(wrap(this->bfgs.hessian)),as<mat>(wrap(hessian0)));
      double edf = trace(Proj);
      NumericVector edf_var = as<NumericVector>(wrap(Smooth::traces(Proj)));
      this->post_process();
      return List::create(_("sp")=wrap(sp),
			  _("coef")=wrap(this->bfgs.coef),
			  _("hessian")=wrap(this->bfgs.hessian),
			  _("edf")=wrap(edf),
			  _("edf_var")=wrap(edf_var)
			  );
    }
    SEXP optim_first() { 
      this->bfgs.reltol = reltol_search;
      double opt_sp = exp(Brent_fmin(log(0.001),log(1000.0),&(pstpm2_first_step<This>),(void *) this,reltol_outer));
      sp[0] = opt_sp;
      this->bfgs.reltol = this->reltol;
      return optim_fixed();
    }
    SEXP optim_multivariate() {
      if (outer_optim==1) return optim_multivariate_NelderMead(); 
      else return optim_multivariate_Nlm(); 
    }
    SEXP optim_multivariate_NelderMead() {
      this->kappa = 10.0;
      NelderMead2 nm;
      nm.reltol = reltol_outer;
      nm.maxit = 500;
      nm.hessianp = false;
      nm.parscale = this->parscale;
      this->bfgs.reltol = reltol_search;
      NumericVector logsp(sp.size());
      for (size_t i=0; i < sp.size(); ++i)
	logsp[i] = log(sp[i]);
      bool satisfied;
      do {
	nm.optim(&pstpm2_multivariate_step<This>, logsp, (void *) this);
	satisfied = true;
	for (size_t i=0; i < sp.size(); ++i)
	  if (logsp[i]< -7.0 || logsp[i] > 7.0) satisfied = false;
	if (!satisfied) this->kappa *= 2.0;
      } while (!satisfied && this->kappa<1.0e5);
      for (int i=0; i < nm.coef.size(); ++i)
	sp[i] = exp(nm.coef[i]);
      this->bfgs.coef = this->init;
      this->bfgs.reltol = this->reltol;
      return optim_fixed();
    }
    SEXP optim_multivariate_Nlm() {
      this->kappa = 10.0;
      Nlm2 nm;
      nm.gradtl = nm.steptl = reltol_outer;
      nm.itnlim = 100;
      // nm.hessianp = false;
      nm.parscale = this->parscale;
      this->bfgs.reltol = reltol_search;
      NumericVector logsp(sp.size());
      for (size_t i=0; i < sp.size(); ++i)
  logsp[i] = log(sp[i]);
      bool satisfied;
      do {
  nm.optim(&pstpm2_multivariate_stepNlm<This>, logsp, (void *) this);
	satisfied = true;
	for (size_t i=0; i < sp.size(); ++i)
	  if (logsp[i]< -7.0 || logsp[i] > 7.0) satisfied = false;
	if (!satisfied) this->kappa *= 2.0;
      } while (!satisfied && this->kappa<1.0e5);
      for (int i=0; i < nm.coef.size(); ++i)
	sp[i] = exp(nm.coef[i]);
      this->bfgs.coef = this->init;
      this->bfgs.reltol = this->reltol;
      return optim_fixed();
    }
    vec sp;
    double alpha, reltol_search, reltol_outer; 
    int criterion, outer_optim;
  };
  // template<class Smooth, class Stpm2>
  // void pstpm2_step_multivariate_nlm(int n, double * logsp, double * f, void *ex) {
  //   *f = pstpm2_step_multivariate<Smooth,Stpm2>(n, logsp, ex);
  // };

  /** 
      Extension of stpm2 and pstpm2 to include gamma shared frailties with a cluster variable
  **/
  template<class Base>
  class GammaSharedFrailty : public Base {
  public:
    typedef std::vector<int> Index;
    typedef std::map<int,Index> IndexMap;
    typedef GammaSharedFrailty<Base> This;
    GammaSharedFrailty(SEXP sexp) : Base(sexp) {
      List list = as<List>(sexp);
      ivec cluster = as<ivec>(list["cluster"]);
      // wragged array indexed by a map of vectors
      for (size_t id=0; id<cluster.size(); ++id) {
	clusters[cluster[id]].push_back(id);
	if (this->event[id]==1)
	  cluster_events[cluster[id]].push_back(id);
      }
    }
    /** 
	Objective function for a Gamma shared frailty
	Assumes that weights == 1.
    **/
    double objective(vec beta) {
      int n = beta.size();
      if (this->bfgs.trace>0) {
	Rprint(beta);
      }
      vec vbeta(beta); // logtheta is the last parameter in beta
      vbeta.resize(n-1);
      double theta = exp(beta[n-1]);
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec h = Base::h(eta,etaD) + this->bhazard;
      vec H = Base::H(eta);
      double constraint = this->kappa/2.0 * (sum(h % h % (h<0)) + sum(H % H % (H<0)));
      vec eps = h*0.0 + 1.0e-16; 
      h = max(h,eps);
      H = max(H,eps);
      vec H0;
      if (this->delayed) {
	vec eta0 = this->X0 * vbeta;
	H0 = Base::H(eta0);
	constraint += this->kappa/2.0 * sum(H0 % H0 % (H0<0));
      }
      double ll = 0.0;
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	uvec index = conv_to<uvec>::from(it->second);
	int mi = sum(this->event(index));
	double sumH = sum(H(index));
	double sumHenter = 0.0;
	ll += sum(log(h(index)) % this->event(index));
	if (this->delayed && cluster_events.count(it->first) > 0) {
	  uvec ind00 = conv_to<uvec>::from(cluster_events[it->first]);
	  sumHenter = sum(H0(ind00));
	}
	 ll += -(1.0/theta+mi)*log(1.0+theta*(sumH)) + 1.0/theta*log(1.0+theta*sumHenter); // Rondeau et al
	// ll -= (1.0/theta+mi)*log(1.0+theta*(sumH - sumHenter)); // conditional (Gutierrez 2002)
	if (mi>0) {
	  for (int k=1; k<=mi; ++k)
	    ll += log(1.0+theta*(mi-k));
	}
      }
      ll -= constraint;
      return -ll;  
    }
    vec gradient(vec beta) {
      int n = beta.size();
      vec gr = zeros<vec>(n);
      vec grconstraint = zeros<vec>(n);
      vec vbeta(beta); // theta is the last parameter in beta
      vbeta.resize(n-1);
      double theta = exp(beta[n-1]);
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec h = Base::h(eta,etaD);
      vec H = Base::H(eta);
      mat gradh = Base::gradh(eta,etaD,this->X,this->XD);
      mat gradH = Base::gradH(eta,this->X);
      vec eps = h*0.0 + 1.0e-16; 
      h = max(h,eps);
      H = max(H,eps);
      vec H0;
      mat gradH0;
      if (this->delayed) {
	vec eta0 = this->X0 * vbeta;
	// vec etaD0 = this->XD0 * vbeta;
	H0 = Base::H(eta0);
	gradH0 = Base::gradH(eta0,this->X0);
      }
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	int mi=0;
	double sumH = 0.0, sumHenter = 0.0;
	vec gradi = zeros<vec>(n-1);
	vec gradHi = zeros<vec>(n-1);
	vec gradH0i = zeros<vec>(n-1);
	vec grconstraint = zeros<vec>(n-1);
	for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	  sumH += H(*j);
	  gradHi += gradH.row(*j).t();
	  gradi += gradh.row(*j).t() / h(*j) * this->event(*j);
	  if (this->event(*j)==1) 
	    mi++;
	  if (this->delayed  && this->ind0[*j]) {
	    int jj = this->map0[*j];
	    sumHenter += H0(jj);
	    gradH0i += gradH0.row(jj).t();
	    grconstraint += this->kappa * (gradH0.row(jj).t() * H0(jj) * (H0(jj)<0));
	  }
	  grconstraint += this->kappa * ((gradh.row(*j).t() * h(*j) * (h(*j)<0)) + (gradH.row(*j).t() * H(*j) * (H(*j)<0)));
	}
	for (int k=0; k<n-1; ++k) {
	  // gr(k) += gradi(k) - (1+mi*theta)*(gradHi(k)-gradH0i(k))/(1+theta*(sumH-sumHenter)) - grconstraint(k); // Gu tierrez
	   gr(k) += gradi(k) - (1+mi*theta)*gradHi(k)/(1+theta*(sumH)) + 1.0/(1+theta*sumHenter)*gradH0i(k) - grconstraint(k); // Rondeau et al
	}
	double lastterm = 0.0;
	if (mi>0) {
	  for (int k=1; k<=mi; ++k)
	    lastterm += theta*(mi-k)/(1.0+theta*(mi-k));
	}
	// gr(n-1) += log(1.0+theta*(sumH-sumHenter))/theta - 
	//   (1.0/theta+mi)*theta*(sumH-sumHenter)/(1.0+theta*(sumH-sumHenter)) + 
	//   lastterm; // Gutierrez 
	gr(n-1) += log(1.0+theta*sumH)/theta - 
	  (1.0+mi*theta)*sumH/(1.0+theta*sumH) + 
	  sumHenter/(1.0+theta*sumHenter) -
	  log(1.0+theta*sumHenter)/theta +
	  lastterm; // Rondeau et al
      }
      if (this->bfgs.trace>1) {
	Rprintf("Calculating fdgradient:\n"); Rprint(this->fdgradient(beta)); Rprintf("(=fdgradient)\n");
      }
      return -gr;  
    }
    bool feasible(vec beta) {
      vec coef(beta);
      coef.resize(this->n-1);
      return Base::feasible(coef);
    }
    void optimWithConstraint(NumericVector init) {
      bool satisfied;
      if (this->bfgs.trace > 0) 
	Rprintf("Starting optimization\n");
      do {
	this->bfgs.optim(&optimfunction<This>, &optimgradient<This>, this->init, (void *) this);
	vec vcoef(&this->bfgs.coef[0],this->n); 
	satisfied = feasible(vcoef % this->parscale);
	if (!satisfied) this->kappa *= 2.0;   
      } while ((!satisfied) && this->kappa < 1.0e3);
    }
    void optimWithConstraintNM(NumericVector init) {
      bool satisfied;
      NelderMead2 nm;
      nm.hessianp = false;
      nm.parscale = this->parscale;
      do {
	nm.optim(&optimfunction<This>, this->init, (void *) this);
	vec vcoef(&this->bfgs.coef[0],this->n); 
	satisfied = feasible(vcoef % this->parscale);
	if (!satisfied) this->kappa *= 2.0;   
      } while ((!satisfied) && this->kappa < 1.0e3);
      nm.hessian = nm.calc_hessian(&optimfunction<This>, (void *) this);
      this->bfgs.coef = nm.coef;
      this->bfgs.hessian = nm.hessian;
    }
    IndexMap clusters, cluster_events;
  };

  /** 
      Extension of stpm2 and pstpm2 to include log-normal shared frailties with a cluster variable
  **/
  template<class Base>
  class NormalSharedFrailty : public Base {
  public:
    typedef std::vector<int> Index;
    typedef std::map<int,Index> IndexMap;
    typedef NormalSharedFrailty<Base> This;
    NormalSharedFrailty(SEXP sexp) : Base(sexp) {
      List list = as<List>(sexp);
      IntegerVector cluster = as<IntegerVector>(list["cluster"]);
      gauss_x = as<vec>(list["gauss_x"]); // node values
      gauss_w = as<vec>(list["gauss_w"]); // probability weights
      // wragged array indexed by a map of vectors
      for (int id=0; id<cluster.size(); ++id) {
	clusters[cluster[id]].push_back(id);
	if (this->event[id]==1) {
	  cluster_events[cluster[id]].push_back(id);
	  cluster_event_indices[cluster[id]].push_back(id);
	}
      }
    }
    /** 
	Objective function for a log-normal shared frailty
	Assumes that weights == 1.
    **/
    double objective(vec beta) {
      int n = beta.size();
      vec vbeta(beta); // nu=log(variance) is the last parameter in beta; exp(nu)=sigma^2
      vbeta.resize(n-1);
      double sigma = exp(0.5*beta[n-1]); // standard deviation
      int K = gauss_x.size(); // number of nodes
      mat lijk(this->N,K);
      double constraint = 0.0;
      vec wstar = gauss_w/sqrt(datum::pi);
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec eta0 = this->X0 * vbeta;
      vec eta1 = this->X1 * vbeta;
      for (int k=0; k<K; ++k) {
	double bi = sqrt(2.0)*sigma*this->gauss_x[k];
	li_constraint lik = Base::li(eta+bi,etaD,eta0+bi,eta1+bi);
	lijk.col(k) = lik.li;
	constraint += lik.constraint;
      }
      double ll = 0.0;
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	uvec index = conv_to<uvec>::from(it->second);
	double Li = dot(exp(sum(lijk.rows(index),0)),wstar);
	ll += log(Li);
      }
      // ll -= constraint;
      return -ll;  
    }
/// Another way for gradients
/// gradient of objective
vec gradient_new(vec beta) {
	int n = beta.size();
	double sigma = exp(0.5*beta[n-1]);
	mat Z(this->N,1,fill::ones);
	mat ZD(this->N,1,fill::zeros);
	mat Z0(this->X0.n_rows,1,fill::ones);
	mat Z1(this->X1.n_rows,1,fill::ones);
	mat Xstar = join_horiz(this->X, Z); 
	mat XDstar = join_horiz(this->XD, ZD); // assumes time-invariant random effects
	mat X0star = join_horiz(this->X0, Z0); 
	mat X1star = join_horiz(this->X1, Z1); 
	int K = gauss_x.size(); // number of nodes
	vec wstar = gauss_w/sqrt(datum::pi);
	int G = clusters.size();
	vec Li(G,fill::zeros);
	mat gradLi(G,n,fill::zeros);
	vec betastar = beta;
	rowvec dconstraint(n,fill::zeros); 
	mat likeli_jk(this->N,K,fill::zeros);
	vec likeli_k(K,fill::zeros);
	cube grad_jk(this->N,n,K,fill::zeros);
	mat grad_kn(K,n,fill::zeros);
	// for K nodes calculation 
	for (int k=0; k<K; ++k) {
		double bi = sqrt(2.0)*sigma*gauss_x[k];
		betastar[betastar.size()-1] = bi;
		vec etastar = Xstar * betastar;
		vec etaDstar = XDstar * betastar;
		vec eta0star = X0star * betastar;
		vec eta1star = X1star * betastar;
		li_constraint lik = Base::li(etastar, etaDstar, eta0star, eta1star);
		gradli_constraint gradlik = Base::gradli(etastar, etaDstar, eta0star, eta1star,Xstar, XDstar, X0star, X1star);
		// adjust the last column of the gradient to account for the variance components:
			// chain rule: d lik/d bi * d bi/d nu where bi = sqrt(2)*exp(nu/2)*x_k
		gradlik.gradli.col(gradlik.gradli.n_cols-1) *= bi*0.5;
		dconstraint += sum(gradlik.constraint,0); // ignores clustering??
		likeli_jk.col(k) = lik.li;
		grad_jk.slice(k) = gradlik.gradli;
	}
	// Iteration in each cluster and get vec Li(g) and mat gradLi(g,n)
	int g = 0;
	for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it, ++g) {
	  uvec index = conv_to<uvec>::from(it->second);
		vec likeli_k = exp(sum(likeli_jk.rows(index),0)).t(); // sum  index components at all nodes for cluster i
		for(int k=0; k<K; ++k){
			grad_kn.row(k) = sum(grad_jk.slice(k).rows(index),0);
		}
		Li(g) = dot(likeli_k,wstar);
		gradLi.row(g) = sum(rmult(grad_kn,wstar%likeli_k),0);
		// gradLi.row(g) = grad_kn.t() * (wstar%likeli_k);
	}
	// Last step to summarize all clusters 
	rowvec grad = sum(rmult(gradLi,1/Li),0); // Rows of matrix is multiplied with elements of vector
	if (this->bfgs.trace>0) {
		Rprintf("fdgradient="); Rprint(this->fdgradient(beta)); 
	}
	vec gr(n,fill::zeros);
	for (int i=0; i<n; i++) gr(i) = grad(i); // - dconstraint(i);
	return -gr;  
}
    vec gradient(vec beta) {
      int n = beta.size();
      double sigma = exp(0.5*beta[n-1]);
      mat Z(this->N,1,fill::ones);
      mat ZD(this->N,1,fill::zeros);
      mat Z0(this->X0.n_rows,1,fill::ones);
      mat Z1(this->X1.n_rows,1,fill::ones);
      mat Xstar = join_horiz(this->X, Z); 
      mat XDstar = join_horiz(this->XD, ZD); // assumes time-invariant random effects
      mat X0star = join_horiz(this->X0, Z0); 
      mat X1star = join_horiz(this->X1, Z1); 
      int K = gauss_x.size(); // number of nodes
      vec wstar = gauss_w/sqrt(datum::pi);
      int G = clusters.size();
      vec Li(G,fill::zeros);
      mat gradLi(G,n,fill::zeros);
      vec betastar = beta;
      rowvec dconstraint(n,fill::zeros); 
      for (int k=0; k<K; ++k) {
	double bi = sqrt(2.0)*sigma*gauss_x[k];
	betastar[betastar.size()-1] = bi;
	vec etastar = Xstar * betastar;
	vec etaDstar = XDstar * betastar;
	vec eta0star = X0star * betastar;
	vec eta1star = X1star * betastar;
	li_constraint lik = Base::li(etastar, etaDstar, eta0star, eta1star);
	gradli_constraint gradlik = Base::gradli(etastar, etaDstar, eta0star, eta1star,
						 Xstar, XDstar, X0star, X1star);
	// adjust the last column of the gradient to account for the variance components:
	// chain rule: d lik/d bi * d bi/d nu where bi = sqrt(2)*exp(nu/2)*x_k
	gradlik.gradli.col(gradlik.gradli.n_cols-1) *= bi*0.5;
	dconstraint += sum(gradlik.constraint,0); // ignores clustering??
	int g = 0;
	for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it, ++g) {
	  uvec index = conv_to<uvec>::from(it->second);
	  double Lgk = exp(sum(lik.li(index)));
	  Li(g) += Lgk*wstar[k];
	  gradLi.row(g) += wstar[k]*Lgk*sum(gradlik.gradli.rows(index),0);
	}
      }
      rowvec grad = sum(rmult(gradLi,1/Li),0);
      if (this->bfgs.trace>0) {
	Rprintf("fdgradient="); Rprint(this->fdgradient(beta)); 
      }
      vec gr(n,fill::zeros);
      for (int i=0; i<n; i++) gr(i) = grad(i); // - dconstraint(i);
      return -gr;  
    }
    bool feasible(vec beta) {
      vec coef(beta);
      coef.resize(beta.size()-1);
      return Base::feasible(coef);
    }
    void optimWithConstraint(NumericVector init) {
      bool satisfied;
      int n = init.size();
      if (this->bfgs.trace > 0) 
	Rprintf("Starting optimization\n");
      do {
	this->bfgs.optim(&optimfunction<This>, &optimgradient<This>, init, (void *) this);
	vec vcoef(&this->bfgs.coef[0],n); 
	satisfied = feasible(vcoef % this->parscale);
	if (!satisfied) this->kappa *= 2.0;   
      } while ((!satisfied) && this->kappa < 1.0e3);
    }
    void optimWithConstraintNM(NumericVector init) {
      bool satisfied;
      NelderMead2 nm;
      nm.hessianp = false;
      nm.parscale = this->parscale;
      do {
	nm.optim(&optimfunction<This>, this->init, (void *) this);
	vec vcoef(&this->bfgs.coef[0],this->n); 
	satisfied = feasible(vcoef % this->parscale);
	if (!satisfied) this->kappa *= 2.0;   
      } while ((!satisfied) && this->kappa < 1.0e3);
      nm.hessian = nm.calc_hessian(&optimfunction<This>, (void *) this);
      this->bfgs.coef = nm.coef;
      this->bfgs.hessian = nm.hessian;
    }
    IndexMap clusters, cluster_events, cluster_event_indices;
    vec gauss_x, gauss_w;
  };



  template<class T>
  SEXP stpm2_model_output_(SEXP args) {
    T model(args);
    List list = as<List>(args);
    std::string return_type = as<std::string>(list["return_type"]);
    vec beta(&model.init[0],model.n);
    if (return_type == "optim") {
      model.pre_process();
      model.optimWithConstraint(model.init);
      model.bfgs.hessian = model.bfgs.calc_hessian(&optimgradient<T>, (void *) &model);
      model.post_process();
      return List::create(_("fail")=model.bfgs.fail, 
			  _("coef")=wrap(model.bfgs.coef),
			  _("hessian")=wrap(model.bfgs.hessian));
    }
    else if (return_type == "objective")
      return wrap(model.objective(beta));
    else if (return_type == "gradient")
      return wrap(model.gradient(beta));
    else if (return_type == "feasible")
      return wrap(model.feasible(beta));
    else {
      REprintf("Unknown return_type.\n");
      return wrap(-1);
    }
  }
  template<class T>
  SEXP pstpm2_model_output_(SEXP args) {
    T model(args);
    List list = as<List>(args);
    std::string return_type = as<std::string>(list["return_type"]);
    vec beta(&model.init[0],model.n);
    if (return_type == "optim_fixed")
      return(model.optim_fixed());
    else if (return_type == "optim_first")
      return(model.optim_first());
    else if (return_type == "optim_multivariate")
      return(model.optim_multivariate());
    else if (return_type == "objective")
      return wrap(model.objective(beta));
    else if (return_type == "objective0")
      return wrap(model.objective0(beta));
    else if (return_type == "gradient")
      return wrap(model.gradient(beta));
    else if (return_type == "gradient0")
      return wrap(model.gradient0(beta));
    else if (return_type == "constraint")
      return wrap(model.feasible(beta));
    else if (return_type == "feasible")
      return wrap(model.feasible(beta));
    else {
      REprintf("Unknown return_type.\n");
      return wrap(-1);
    }
  }
  RcppExport SEXP model_output(SEXP args) {
    List list = as<List>(args);
    std::string link = as<std::string>(list["link"]);
    std::string type = as<std::string>(list["type"]);
    if (type=="stpm2") {
      if (link == "PH")
	return stpm2_model_output_<Stpm2<LogLogLink> >(args);
      else if (link == "PO")
	return stpm2_model_output_<Stpm2<LogitLink> >(args);
      else if (link == "AH")
	return stpm2_model_output_<Stpm2<LogLink> >(args);
      else if (link == "probit")
	return stpm2_model_output_<Stpm2<ProbitLink> >(args);
      else {
	REprintf("Unknown link function.\n");
	return wrap(-1);
      }
    }
    else if (type=="pstpm2") {
      if (link == "PH")
	return pstpm2_model_output_<Pstpm2<Stpm2<LogLogLink>,SmoothLogH> >(args);
      else if (link == "PO")
	return pstpm2_model_output_<Pstpm2<Stpm2<LogitLink>,SmoothLogH> >(args);
      else if (link == "AH")
	return pstpm2_model_output_<Pstpm2<Stpm2<LogLink>,SmoothLogH> >(args);
      else if (link == "probit")
	return pstpm2_model_output_<Pstpm2<Stpm2<ProbitLink>,SmoothLogH> >(args);
      else {
	REprintf("Unknown link function.\n");
	return wrap(-1);
      }
    }
    else if (type=="stpm2_gamma_frailty") {
      if (link == "PH")
	return stpm2_model_output_<GammaSharedFrailty<Stpm2<LogLogLink> > >(args);
      else if (link == "PO")
	return stpm2_model_output_<GammaSharedFrailty<Stpm2<LogitLink> > >(args);
      else if (link == "AH")
	return stpm2_model_output_<GammaSharedFrailty<Stpm2<LogLink> > >(args);
      else if (link == "probit")
	return stpm2_model_output_<GammaSharedFrailty<Stpm2<ProbitLink> > >(args);
      else {
	REprintf("Unknown link function.\n");
	return wrap(-1);
      }
    }
    else if (type=="pstpm2_gamma_frailty") {
      if (link == "PH")
	return pstpm2_model_output_<Pstpm2<GammaSharedFrailty<Stpm2<LogLogLink> >,SmoothLogH> >(args);
      else if (link == "PO")
	return pstpm2_model_output_<Pstpm2<GammaSharedFrailty<Stpm2<LogitLink> >,SmoothLogH> >(args);
      else if (link == "AH")
	return pstpm2_model_output_<Pstpm2<GammaSharedFrailty<Stpm2<LogLink> >,SmoothLogH> >(args);
      else if (link == "probit")
	return pstpm2_model_output_<Pstpm2<GammaSharedFrailty<Stpm2<ProbitLink> >,SmoothLogH> >(args);
      else {
	REprintf("Unknown link function.\n");
	return wrap(-1);
      }
    }
    else if (type=="stpm2_normal_frailty") {
      if (link == "PH")
	return stpm2_model_output_<NormalSharedFrailty<Stpm2<LogLogLink> > >(args);
      else if (link == "PO")
	return stpm2_model_output_<NormalSharedFrailty<Stpm2<LogitLink> > >(args);
      else if (link == "AH")
	return stpm2_model_output_<NormalSharedFrailty<Stpm2<LogLink> > >(args);
      else if (link == "probit")
	return stpm2_model_output_<NormalSharedFrailty<Stpm2<ProbitLink> > >(args);
      else {
	REprintf("Unknown link function.\n");
	return wrap(-1);
      }
    }
    else if (type=="pstpm2_normal_frailty") {
      if (link == "PH")
	return pstpm2_model_output_<Pstpm2<NormalSharedFrailty<Stpm2<LogLogLink> >,SmoothLogH> >(args);
      else if (link == "PO")
	return pstpm2_model_output_<Pstpm2<NormalSharedFrailty<Stpm2<LogitLink> >,SmoothLogH> >(args);
      else if (link == "AH")
	return pstpm2_model_output_<Pstpm2<NormalSharedFrailty<Stpm2<LogLink> >,SmoothLogH> >(args);
      else if (link == "probit")
	return pstpm2_model_output_<Pstpm2<NormalSharedFrailty<Stpm2<ProbitLink> >,SmoothLogH> >(args);
      else {
	REprintf("Unknown link function.\n");
	return wrap(-1);
      }
    }
    else {
      REprintf("Unknown model type.\n");
      return wrap(-1);
    }
  }


  // Proof of concept for a Weibull cure model
  
  struct CureModel {
    int n0, n1, n2;
    mat Xshape, Xscale, Xcure;
    vec time, status;
  };

  double fminfn_cureModel(int n, double * beta, void *ex) {
    double ll = 0.0;
    CureModel * data = (CureModel *) ex;
    vec vbeta(&beta[0],n);
    vec shape = exp(data->Xshape * vbeta(span(0,data->n0-1)));
    vec scale = exp(data->Xscale * vbeta(span(data->n0,data->n1-1)));
    vec cure = 1.0/(1.0+exp(-data->Xcure * vbeta(span(data->n1,data->n2-1))));
    for (size_t i=0; i < data->time.size(); ++i) {
      ll += data->status(i)==1.0 ?
	log(1.0-cure(i)) + R::dweibull(data->time(i),shape(i),scale(i),1) :
	log(cure(i)+(1.0-cure(i)) * R::pweibull(data->time(i),shape(i),scale(i),0,0));
    }
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
    return -ll;
  }

  RcppExport SEXP fitCureModel(SEXP stime, SEXP sstatus, SEXP sXshape,
			  SEXP sXscale, SEXP sXcure, SEXP sbeta) {
    mat Xshape = as<mat>(sXshape);
    mat Xscale = as<mat>(sXscale);
    mat Xcure = as<mat>(sXcure);
    vec time = as<vec>(stime);
    vec status = as<vec>(sstatus);
    NumericVector init = as<NumericVector>(sbeta);
    int n0=Xshape.n_cols;
    int n1=n0+Xscale.n_cols;
    int n2=n1+Xcure.n_cols;
    CureModel data = {n0,n1,n2,Xshape,Xscale,Xcure,time,status};
    
    NelderMead nm;
    nm.reltol = 1.0e-8;
    nm.maxit = 1000;
    nm.optim(& fminfn_cureModel, init, (void *) &data);

    for (int i = 0; i<nm.coef.size(); ++i)
      init[i] = nm.coef[i]; // clone
    
    return List::create(_("Fmin")=nm.Fmin, 
			_("coef")=wrap(init),
			_("fail")=nm.fail,
			_("hessian")=wrap(nm.hessian));
  }

  // R CMD INSTALL ~/src/R/microsimulation
  // R -q -e "require(microsimulation); .Call('test_nmmin',1:2,PACKAGE='microsimulation')"

} // anonymous namespace
