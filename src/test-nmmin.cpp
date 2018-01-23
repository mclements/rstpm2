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
  // constants
  const double log2pi = std::log(2.0 * M_PI); // 
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
  void Rprint(uvec const & v) {
    for (size_t i=0; i<v.size(); ++i) 
      Rprintf("%i ", v(i));
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
  void Rprint(cube const & m) {
    for (size_t k=0; k<m.n_slices; ++k) {
      Rprintf("[");
      for (size_t i=0; i<m.n_rows; ++i) {
	for (size_t j=0; j<m.n_cols; ++j) 
	  Rprintf("%f ", m(i,j,k));
	Rprintf("\n");
      }
      Rprintf("]\n");
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
  // http://gallery.rcpp.org/articles/dmvnorm_arma/
  arma::vec dmvnrm_arma(arma::mat x,  
			arma::rowvec mean,  
			arma::mat sigma, 
			bool logd = false) { 
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    for (int i=0; i < n; i++) {
      arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
      out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
    }  
    if (logd == false) {
      out = exp(out);
    }
    return(out);
  }
  double dmvnrm_arma(arma::vec x,  
		     arma::vec mean,  
		     arma::mat sigma, 
		     bool logd = false) { 
    int xdim = x.n_rows;
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    arma::vec z = rooti * (x - mean) ;    
    double out = constants - 0.5 * arma::sum(z%z) + rootisum;     
    if (logd == false) {
      out = exp(out);
    }
    return(out);
  }
  
  class Link {
  public:
    virtual vec link(vec S) = 0;
    virtual vec ilink(vec x) = 0;
    virtual vec h(vec eta, vec etaD) = 0;
    virtual vec H(vec eta) = 0;
    virtual mat gradH(vec eta, mat X) = 0;
    virtual mat gradh(vec eta, vec etaD, mat X, mat XD) = 0;
    virtual ~Link() { }
  };
  // Various link objects
  class LogLogLink : public Link {
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
    LogLogLink(SEXP sexp) {}
  };
  class ArandaOrdazLink : public Link {
  public:
    vec link(vec S) {
      if (thetaAO==0)
	return log(-log(S));
      else return log((exp(log(S)*(-thetaAO))-1)/thetaAO);
    }
    vec ilink(vec eta) {
      if (thetaAO==0)
	return exp(-exp(eta));
      else return exp(-log(thetaAO*exp(eta)+1)/thetaAO);
    }
    vec h(vec eta, vec etaD) {
      if (thetaAO==0)
	return etaD % exp(eta);
      else return exp(eta) % etaD/(thetaAO*exp(eta)+1);
    }
    vec H(vec eta) {
      if (thetaAO==0)
	return exp(eta);
      else return log(thetaAO*exp(eta)+1)/thetaAO;
    }
    mat gradH(vec eta, mat X) {
      if (thetaAO==0)
	return rmult(X,exp(eta));
      else return rmult(X,exp(eta)/(1+thetaAO*exp(eta)));
    }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      if (thetaAO==0)
	return rmult(XD, exp(eta)) + rmult(X, etaD % exp(eta));
      else {
	vec denom = (thetaAO*exp(eta)+1) % (thetaAO*exp(eta)+1);
	return rmult(XD,(thetaAO*exp(2*eta)+exp(eta))/denom)+rmult(X,exp(eta) % etaD/denom);
      }
    }
    double thetaAO;
    ArandaOrdazLink(SEXP sexp) { 
      List args = as<List>(sexp);
      thetaAO = as<double>(args["thetaAO"]);
    }
  };
  // Useful relationship: d/dx expit(x)=expit(x)*expit(-x) 
  class LogitLink : public Link {
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
    LogitLink(SEXP sexp) {}
  };
  class ProbitLink : public Link {
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
    // incomplete implementation
    cube hessianh(vec beta, mat X, mat XD) { 
      cube c(beta.size(), beta.size(), X.n_rows, fill::zeros);
      return c;
    }
    // incomplete implementation
    cube hessianH(vec beta, mat X) { 
      cube c(beta.size(), beta.size(), X.n_rows, fill::zeros);
      return c;
    }
    ProbitLink(SEXP sexp) {}
  };
  class LogLink : public Link {
  public:
    vec link(vec S) { return -log(S); }
    vec ilink(vec x) { return exp(-x); }
    vec h(vec eta, vec etaD) { return etaD; }
    vec H(vec eta) { return eta; }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { return XD; }
    mat gradH(vec eta, mat X) { return X;  }
    cube hessianh(vec beta, mat X, mat XD) { 
      cube c(beta.size(), beta.size(), X.n_rows, fill::zeros);
      return c;
    }
    cube hessianH(vec beta, mat X) { 
      cube c(beta.size(), beta.size(), X.n_rows, fill::zeros);
      return c;
    }
    LogLink(SEXP sexp) {}
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
    // R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
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
    // if (obj->bfgs.trace>2) {
    //   Rprintf("fdgradient="); Rprint(obj->fdgradient(coef % obj->parscale));
    // }
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

  /** @brief Macro to define the optimisation functions. Defined as a macro
      because the C API needs to pass a 'This' template parameter
      which varies for the different classes (alternatively, this
      could be defined using the curiously recurring template
      pattern).
  **/
#define OPTIM_FUNCTIONS							\
    void optimWithConstraintBFGS(NumericVector init) {			\
      bool satisfied;							\
      this->bfgs.coef = init;						\
      if (this->robust_initial) {					\
	optimNoHessianNM(this->bfgs.coef,50);				\
      }									\
      this->kappa = this->kappa_init;					\
      do {								\
	this->bfgs.optim(&optimfunction<This>, &optimgradient<This>, this->bfgs.coef, (void *) this); \
	vec vcoef(&this->bfgs.coef[0],this->n);				\
	satisfied = this->feasible(vcoef % this->parscale);		\
	if (!satisfied) this->kappa *= 2.0;				\
      } while ((!satisfied) && this->kappa < 1.0e3);			\
      if (this->bfgs.trace > 0 && this->kappa>1.0) Rprintf("kappa=%f\n",this->kappa); \
      this->bfgs.hessian = this->bfgs.calc_hessian(&optimgradient<This>, (void *) this); \
    }									\
    void optimNoHessianNM(NumericVector init, int maxit = 50) {		\
      NelderMead2 nm;							\
      nm.hessianp = false;						\
      nm.parscale = this->parscale;					\
      nm.maxit = maxit;							\
      nm.optim(&optimfunction<This>, init, (void *) this);		\
      this->bfgs.coef = nm.coef;					\
    }									\
    void optimWithConstraintNM(NumericVector init) {			\
      bool satisfied;							\
      NelderMead2 nm;							\
      nm.hessianp = false;						\
      nm.parscale = this->parscale;					\
      this->kappa = this->kappa_init;					\
      do {								\
	nm.optim(&optimfunction<This>, init, (void *) this);		\
	vec vcoef(&nm.coef[0],this->n);					\
	satisfied = this->feasible(vcoef % this->parscale);		\
	if (!satisfied) this->kappa *= 2.0;				\
      } while ((!satisfied) && this->kappa < this->maxkappa);			\
      nm.hessian = nm.calc_hessian(&optimfunction<This>, (void *) this); \
      this->bfgs.coef = nm.coef;					\
      this->bfgs.hessian = nm.hessian;					\
    }									\
    void optimWithConstraintNlm(NumericVector init) {			\
      bool satisfied;							\
      Nlm2 nm;								\
      nm.gradtl = nm.steptl = this->reltol;				\
      nm.parscale = this->parscale;					\
      this->kappa = this->kappa_init;					\
      do {								\
	nm.optim(&optimfunction_nlm<This>, init, (void *) this);	\
	vec vcoef(&nm.coef[0],this->n);					\
	satisfied = this->feasible(vcoef % this->parscale);		\
	if (!satisfied) this->kappa *= 2.0;				\
      } while ((!satisfied) && this->kappa < 1.0e3);			\
      if (this->bfgs.trace > 0 && this->kappa>1.0) Rprintf("kappa=%f\n",this->kappa); \
      nm.hessian = nm.calc_hessian(&optimfunction_nlm<This>, (void *) this); \
      this->bfgs.coef = nm.coef;					\
      this->bfgs.hessian = nm.hessian;					\
    }									\
    void optimWithConstraint(NumericVector init) {			\
      if (this->bfgs.trace > 0) Rprintf("Starting optimization\n");	\
      if (this->optimiser == "NelderMead")				\
	optimWithConstraintNM(init);					\
      else if (this->optimiser == "Nlm")				\
	optimWithConstraintNlm(init);					\
      else								\
	optimWithConstraintBFGS(init);					\
    }
    
  // struct for data
  struct BaseData {
    mat X, XD, X0, X1; 
    vec bhazard,wt,wt0,event,time;
    vec map0;
    uvec ind0, which0;
  };
  // Main class for Stpm2 (fully parametric)
  class Stpm2 : public BaseData {
  public:
    typedef Stpm2 This;
    Stpm2(SEXP sexp) : bfgs() {
      List list = as<List>(sexp);
      std::string linkType = as<std::string>(list["link"]);
      if (linkType=="PH")
	link=new LogLogLink(sexp);
      else if (linkType=="PO")
	link=new LogitLink(sexp);
      else if (linkType=="probit")
	link=new ProbitLink(sexp);
      else if (linkType=="AH")
	link=new LogLink(sexp);
      else if (linkType=="AO")
	link=new ArandaOrdazLink(sexp);
      else {
	REprintf("No matching link.type");
	return;
      }
      bfgs.coef = init = as<NumericVector>(list["init"]);
      X = as<mat>(list["X"]); 
      XD = as<mat>(list["XD"]); 
      bhazard = as<vec>(list["bhazard"]);
      wt = as<vec>(list["wt"]);
      event = as<vec>(list["event"]);
      time = as<vec>(list["time"]);
      delayed = as<bool>(list["delayed"]);
      interval = as<bool>(list["interval"]);
      n = nbeta = init.size(); // number of parameters
      N = X.n_rows; // number of observations
      X1 = as<mat>(list["X1"]);
      X0 = as<mat>(list["X0"]);
      wt0.set_size(N); wt0.fill(0.0);
      if (delayed) {
	wt0 = as<vec>(list["wt0"]);
      }
      // if (interval) {
      // 	X1 = as<mat>(list["X1"]);
      // }
      map0 = as<vec>(list["map0"]); // length N map for X->X0
      full_which0 = as<vec>(list["which0"]); // length N map for X0->X
      ind0 = as<uvec>(list["ind0"]); // length N boolean for X0
      ttype = as<vec>(list["ttype"]); 
      kappa_init = kappa = as<double>(list["kappa"]);
      maxkappa = as<double>(list["maxkappa"]);
      optimiser = as<std::string>(list["optimiser"]);
      bfgs.trace = as<int>(list["trace"]);
      reltol = bfgs.reltol = as<double>(list["reltol"]);
      robust_initial = as<bool>(list["robust_initial"]);
      bfgs.hessianp = false;
      bfgs.parscale = parscale = as<vec>(list["parscale"]);
      which0 = removeNaN(full_which0);
    }
    // log-likelihood components and constraints
    // Note: the li and gradli components depend on other variables (bhazard, wt, wt0, wt, event);
    //       calculations should *not* be done on subsets
    virtual ~Stpm2() {
      delete link;
    }
    li_constraint li_right_censored(vec eta, vec etaD) {
      vec h = link->h(eta,etaD) + bhazard;
      vec H = link->H(eta);
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
      vec H0 = link->H(eta0);
      double constraint = kappa/2.0 * sum(H0 % H0 % (H0<0));
      vec eps = H0*0.0 + 1.0e-16; 
      vec li = wt0 % max(H0,eps);
      li_constraint out = {li, constraint};
      return out;
    }
    li_constraint li_interval(vec eta, vec etaD, vec eta1) {
      vec H = link->H(eta);
      vec H1 = link->H(eta1);
      vec h = link->h(eta, etaD) + bhazard;
      double constraint = 0.0;
      vec eps = H*0.0 + 1e-16;
      H = max(H,eps);
      H1 = max(H1,eps);
      h = max(h,eps);
      vec li(N,fill::zeros);
      uvec index;
      index = find(ttype == 0); // right censored
      if (any(index)) {
	li(index) -= wt(index) % H(index);
	constraint += kappa/2.0 * sum(H(index) % H(index) % (H(index)<0));
      }
      index = find(ttype == 1); // exact
      if (any(index)) {
	li(index) += wt(index) % (log(h(index)) - H(index));
	constraint += kappa/2.0 * (sum(H(index) % H(index) % (H(index)<0))+
				   sum(h(index) % h(index) % (h(index)<0)));
      }
      index = find(ttype == 2); // left censored
      if (any(index)) {
	li(index) += wt(index) % log(1-exp(-H(index)));
	constraint += kappa/2.0 * sum(H(index) % H(index) % (H(index)<0));
      }
      index = find(ttype == 3); // interval censored
      if (any(index)) {
	li(index) += wt(index) % log(exp(-H(index)) - exp(-H1(index)));
	constraint += kappa/2.0 * (sum(H(index) % H(index) % (H(index)<0))+
				   sum(H1(index) % H1(index) % (H1(index)<0)));
      }
      li_constraint out = {li, constraint};
      return out;
    }
    li_constraint li(vec eta, vec etaD, vec eta0, vec eta1) {
      if (interval) {
	return li_interval(eta, etaD, eta1);
      }
      else {
	li_constraint s = li_right_censored(eta, etaD);
	if (delayed && !eta0.empty()) {
	  li_constraint s0 = li_left_truncated(eta0);
	  s.constraint += s0.constraint;
	  if (bfgs.trace > 0) {
	    Rprint(which0);
	  }
	  s.li(which0) += s0.li;
	}
	return s;
      }
    }
    vec getli(vec beta) {
      vec vbeta = beta;
      vbeta.resize(nbeta);
      li_constraint lic = li(X*vbeta, XD*vbeta, X0*vbeta, X1*vbeta);
      return lic.li;
    }
    mat getgradli(vec beta) {
      vec vbeta = beta;
      vbeta.resize(nbeta);
      gradli_constraint gradlic = gradli(X*vbeta, XD*vbeta, X0*vbeta, X1*vbeta, X, XD, X0, X1);
      return gradlic.gradli;
    }
    // negative log-likelihood
    double objective(vec beta) {
      vec vbeta = beta;
      vbeta.resize(nbeta);
      li_constraint s = li(X * vbeta, XD * vbeta, X0 * vbeta, X1 * vbeta);
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
    double f(vec beta, int i, double delta) {
      double betai = beta(i);
      beta(i) += delta;
      double y = objective(beta);
      beta(i) = betai;
      return y;
    }
    double f2(vec beta, int i, int j, double epsi, double epsj) {
      double betai = beta(i), betaj = beta(j);
      beta(i) += epsi;
      beta(j) += epsj;
      double y = objective(beta);
      beta(i) = betai;
      beta(j) = betaj;
      return y;
    }
    mat fdhessian(vec beta) {
      int n = beta.size();
      mat hessian(n,n);
      double eps0 = pow(datum::eps,1.0/3.0);
      double epsi, epsj;
      for (int i=0; i<n; ++i) {
	for (int j=i; j<n; ++j) {
	  if (i==j) {
	    epsi = eps0*(1+std::abs(beta(i)));
	    hessian(i,i) = (-f(beta,i,2.0*epsi)+16.0*f(beta,i,epsi)-30.0*f(beta,i,0.0)+16.0*f(beta,i,-epsi)-
		 f(beta,i,-2.0*epsi))/12.0/epsi/epsi;
	  } else {
	    epsi = eps0*(1+std::abs(beta(i)));
	    epsj = eps0*(1+std::abs(beta(j)));
	    hessian(i,j) = hessian(j,i) = (f2(beta,i,j,epsi,epsj)-f2(beta,i,j,epsi,-epsj)-f2(beta,i,j,-epsi,epsj)+f2(beta,i,j,-epsi,-epsj))/4.0/epsi/epsj; 
	  }
	}
      }
      return hessian;
    }
    // log-likelihood gradient components
    // Note: gradH and gradh are currently calculated at eta and etaD and may not satisfy the constraints, 
    //       while H and h may be transformed for the constraints 
    gradli_constraint gradli_right_censored(vec eta, vec etaD, mat X, mat XD) {
      vec h = link->h(eta,etaD) + bhazard;
      vec H = link->H(eta);
      vec eps = h*0.0 + 1.0e-16;
      mat gradH = link->gradH(eta,X);
      mat gradh = link->gradh(eta,etaD,X,XD);
      mat Xconstraint = kappa * (rmult(gradh, h % (h<0)) +
				 rmult(gradH, H % (H<0)));
      // h = max(h,eps);
      mat Xgrad = -rmult(gradH, wt) + rmult(gradh, event / h % wt);
      gradli_constraint out = {Xgrad, Xconstraint};
      return out;
    }
    gradli_constraint gradli_left_truncated(vec eta0, mat X0) {
      mat gradH0 = link->gradH(eta0, X0); 
      vec H0 = link->H(eta0); 
      vec eps = H0*0.0 + 1.0e-16;
      mat Xconstraint = kappa * rmult(gradH0, H0 % (H0<0));
      mat Xgrad = rmult(gradH0, wt0);
      gradli_constraint out = {Xgrad, Xconstraint};
      return out;
    }
    gradli_constraint gradli_interval_censored(vec eta, vec etaD, vec eta1, 
						       mat X, mat XD, mat X1) {
      vec H = link->H(eta);
      vec h = link->h(eta,etaD);
      vec H1 = link->H(eta1);
      mat gradH = link->gradH(eta,X);
      mat gradH1 = link->gradH(eta1,X1);
      mat gradh = link->gradh(eta, etaD, X, XD);
      vec eps = H*0.0 + 1e-16;
      // mat Xconstraint = kappa * (rmult(gradH, H % (H<eps))+
      // 				 rmult(gradh, h % (h<eps))+
      // 				 rmult(gradH1, H1 % (H1<eps)));
      mat Xconstraint = X * 0.0;
      H = max(H,eps);
      H1 = max(H1,eps);
      h = max(h,eps);
      mat li(N,n,fill::zeros);
      uvec index;
      index = find(ttype == 0); // right censored
      if (any(index)) {
	li.rows(index) -= rmult(gradH.rows(index),wt(index));
	Xconstraint.rows(index) += kappa*(rmult(gradH.rows(index), H(index) % (H(index)<eps(index))));
      }
      index = find(ttype == 1); // exact
      if (any(index)) {
	li.rows(index) += rmult(gradh.rows(index),wt(index) / h(index)) - rmult(gradH.rows(index),wt(index));
	Xconstraint.rows(index) += kappa*(rmult(gradH.rows(index), H(index) % (H(index)<eps(index))) +
					  rmult(gradh.rows(index), h(index) % (h(index)<eps(index))));
      }
      index = find(ttype == 2); // left censored
      if (any(index)) {
	li.rows(index) += rmult(-gradH.rows(index),-exp(-H(index)) / (1-exp(-H(index))) % wt(index));
	Xconstraint.rows(index) += kappa*(rmult(gradH.rows(index), H(index) % (H(index)<eps(index))));
      }
      index = find(ttype == 3); // interval censored
      if (any(index)) {
	vec V = wt(index) / (exp(-H(index)) - exp(-H1(index)));
	li.rows(index) += rmult(gradH1.rows(index),V % exp(-H1(index))) - 
	  rmult(gradH.rows(index),V % exp(-H(index)));
	Xconstraint.rows(index) += kappa*(rmult(gradH.rows(index), H(index) % (H(index)<eps(index))) +
					  rmult(gradH1.rows(index), H1(index) % (H1(index)<eps(index))));
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
	gr[i] = vgr[i];// - dconstraint[i];
      }
      return -gr;
    }
    bool feasible(vec beta) {
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec h = link->h(eta, etaD) + bhazard;
      vec H = link->H(eta);
      bool condition = all((h>0) % (H>0));
      if (delayed) {
	vec eta0 = X0 * beta;
	// vec etaD0 = XD0 * beta;
	vec H0 = link->H(eta0);
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
    OPTIM_FUNCTIONS;
    // NumericMatrix calcHessian(NumericVector nvbeta) {
    //   // not defined for interval-censored or probit link
    //   vec beta(&nvbeta[0],n);
    //   vec eta = X * beta;
    //   vec etaD = XD * beta;
    //   vec h = link->h(eta,etaD);
    //   mat gradh = link->gradh(eta,etaD,X,XD);
    //   cube hessianh = link->hessianh(beta,X,XD);
    //   cube hessianH = link->hessianH(beta,X);
    //   cube hessianH0 = link->hessianH(beta,X0);
    //   mat hessian = -arma::sum(hessianH,2);
    //   if (delayed)
    // 	hessian += arma::sum(hessianH0,2);
    //   for (int i=0; i<N; ++i) 
    // 	hessian += event(i)/h(i)*hessianh.slice(i) -
    // 	  event(i)/h(i)/h(i)*gradh.row(i).t()*gradh.row(i);
    //   return wrap(-hessian);
    // }
    uvec map0f(uvec index) {
      return removeNaN(this->map0(index));
    }
    uvec which0f(uvec index) {
      vec which = this->full_which0(index);
      return find(which == which);
    }
    uvec removeNaN(vec x) {
      vec newx = x;
      uvec index = find(newx == newx); // NaN != NaN
      newx = newx(index);
      return conv_to<uvec>::from(newx);
    }
    SEXP return_modes() { return wrap(-1); } // used in NormalSharedFrailties
    SEXP return_variances() { return wrap(-1); } // used in NormalSharedFrailties
    NumericVector init;
    vec parscale, ttype, full_which0;
    double kappa_init, kappa, maxkappa, reltol;
    bool delayed, interval, robust_initial;
    int n, N, nbeta;
    BFGS2 bfgs;
    std::string optimiser;
    Link * link;
  };

  // penalised smoothers 
  // _Not_ sub-classes of Pstpm2, hence the longer function signatures
  class SmoothLogH {
  public:
    struct Smoother {
      int first_para, last_para;
      mat S;
      bool first;
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
	  as<mat>(lsmoothiS[0]),
	  true
	};
	smooth.push_back(smoothi);
	if (lsmoothiS.size()==2) {
	  Smoother smoothi2 =  {
	    as<int>(lsmoothi["first.para"]) - 1, 
	    as<int>(lsmoothi["last.para"]) - 1, 
	    as<mat>(lsmoothiS[1]),
	    false
	  };
	  smooth.push_back(smoothi2);
	}
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
	// if (!smoothi.first) values[i]=0.0; // edf for second smoother matrix: ignore
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
  template<class Stpm2Type = Stpm2, class Smooth = SmoothLogH>
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
    OPTIM_FUNCTIONS;
    double calc_edf(NumericMatrix hessian0) {
      double max_edf = this->bfgs.hessian.ncol();
      mat m;
      bool flag = arma::solve(m, as<mat>(this->bfgs.hessian), as<mat>(hessian0));
      double edf = flag ? arma::trace(m) : (max_edf*2.0);
      if (edf<0.0) edf = max_edf*2.0;
      return edf;
    }
    double first_step(double logsp) {
      sp[0] = exp(logsp);
      this->pre_process();
      this->optimWithConstraint(this->init);
      // this->bfgs.hessian = this->bfgs.calc_hessian(&optimgradient<This>, (void *) this);
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
      double edf = calc_edf(hessian0);
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
      double edf = calc_edf(hessian0);
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
			  _("edf_var")=wrap(edf_var),
			  _("kappa")=wrap(this->kappa)
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
	  if (logsp[i]< -9.0 || logsp[i] > 9.0) satisfied = false;
	if (!satisfied) this->kappa *= 2.0;
      } while (!satisfied && this->kappa<1.0e5);
      for (int i=0; i < nm.coef.size(); ++i)
	sp[i] = exp(nm.coef[i]);
      this->bfgs.coef = this->init;
      this->bfgs.reltol = this->reltol;
      return optim_fixed();
    }
    SEXP optim_multivariate_Nlm() {
      this->kappa = this->kappa_init;
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
	  if (logsp[i]< -9.0 || logsp[i] > 9.0) satisfied = false;
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
      recurrent = as<bool>(list["recurrent"]);
      this->nbeta = this->n - 1;
      // wragged array indexed by a map of vectors
      for (size_t id=0; id<cluster.size(); ++id) {
	clusters[cluster[id]].push_back(id);
	if (this->delayed && this->ind0[id])
	  cluster_delayed[cluster[id]].push_back(id);
      }
    }
    /** 
	Objective function for a Gamma shared frailty
	Assumes that weights == 1.
    **/
    double objective(vec beta) {
      vec vbeta = beta;
      vbeta.resize(this->nbeta);
      // copy across to this->bfgs.coef?
      li_constraint s = li(this->X * vbeta, this->XD * vbeta, this->X0 * vbeta, this->X1 * vbeta);
      return -sum(s.li) + s.constraint;
    }
    vec getli(vec beta) {
      vec vbeta = beta;
      vbeta.resize(this->nbeta);
      li_constraint lic = li(this->X*vbeta, this->XD*vbeta, this->X0*vbeta, this->X1*vbeta);
      return lic.li;
    }
    mat getgradli(vec beta) {
      vec vbeta = beta;
      vbeta.resize(this->nbeta);
      gradli_constraint gradlic = gradli(this->X*vbeta, this->XD*vbeta, this->X0*vbeta, this->X1*vbeta, this->X, this->XD, this->X0, this->X1);
      return gradlic.gradli;
    }
    li_constraint li(vec eta, vec etaD, vec eta0, vec eta1) {
      vec ll(clusters.size(), fill::zeros);
      vec beta = as<vec>(this->bfgs.coef);
      int n = beta.size();
      vec vbeta(beta); // logtheta is the last parameter in beta
      vbeta.resize(this->nbeta);
      double theta = exp(beta[n-1]);
      eta = this->X * vbeta;
      etaD = this->XD * vbeta;
      vec h = this->link->h(eta,etaD) + this->bhazard;
      vec H = this->link->H(eta);
      double constraint = this->kappa/2.0 * (sum(h % h % (h<0)) + sum(H % H % (H<0)));
      vec eps = h*0.0 + 1.0e-16; 
      h = max(h,eps);
      H = max(H,eps);
      vec H0;
      if (this->delayed) {
	eta0 = this->X0 * vbeta;
	H0 = this->link->H(eta0);
	constraint += this->kappa/2.0 * sum(H0 % H0 % (H0<0));
      }
      int i=0;
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it, ++i) {
	uvec index = conv_to<uvec>::from(it->second);
	int mi = sum(this->event(index));
	double sumH = sum(H(index));
	double sumHenter = 0.0;
	ll(i) += sum(log(h(index)) % this->event(index)); // shared
	if (this->delayed && !cluster_delayed[it->first].empty()) {
	  uvec ind00 = conv_to<uvec>::from(cluster_delayed[it->first]);
	  sumHenter = sum(H0(Base::map0f(ind00)));
	}
	if (recurrent) { // for calendar timescale
		ll(i) += -(1.0/theta+mi)*log(1.0+theta*(sumH - sumHenter)); // Rondeau et al (2005)
	} else { // for clustered data or gap timescale
		ll(i) += -(1.0/theta+mi)*log(1.0+theta*sumH) + 1.0/theta*log(1.0+theta*sumHenter); // Rondeau et al
	}
	// Note: log(gamma(1/theta+mi)/gamma(1/theta)*theta^mi) = sum_k=1^mi (log(1+theta(mi-k))) (mi>0); 0 (otherwise)
	if (mi>0) { 
		for (int k=1; k<=mi; ++k)
			ll(i) += log(1.0+theta*(mi-k));
	} 
      }
      li_constraint out;
      out.constraint=constraint;
      out.li = ll;
      return out;
    }
    vec gradient(vec beta) {
      vec vbeta = beta;
      vbeta.resize(this->nbeta);
      gradli_constraint gc = gradli(this->X*vbeta, this->XD*vbeta, this->X0*vbeta, this->X1*vbeta, this->X, this->XD, this->X0, this->X1);
      rowvec dconstraint = sum(gc.constraint,0);
      rowvec vgr = sum(gc.gradli,0);
      vec gr(beta.size());
      for (size_t i = 0; i<beta.size(); ++i) {
	gr[i] = vgr[i];// - dconstraint[i];
      }
      return -gr;
    }
    gradli_constraint gradli(vec eta, vec etaD, vec eta0, vec eta1, mat X, mat XD, mat X0, mat X1) { 
      vec beta = as<vec>(this->bfgs.coef);
      int n = beta.size();
      mat gr = zeros<mat>(clusters.size(), n);
      mat grconstraint = zeros<mat>(clusters.size(), n);
      vec vbeta(beta); // theta is the last parameter in beta
      vbeta.resize(n-1);
      double theta = exp(beta[n-1]);
      // eta = this->X * vbeta;
      // etaD = this->XD * vbeta;
      vec h = this->link->h(eta,etaD);
      vec H = this->link->H(eta);
      mat gradh = this->link->gradh(eta,etaD,this->X,this->XD);
      mat gradH = this->link->gradH(eta,this->X);
      vec eps = h*0.0 + 1.0e-16; 
      h = max(h,eps);
      H = max(H,eps);
      vec H0;
      mat gradH0;
      if (this->delayed) {
	// vec eta0 = this->X0 * vbeta;
	// vec etaD0 = this->XD0 * vbeta;
	H0 = this->link->H(eta0);
	gradH0 = this->link->gradH(eta0,this->X0);
      }
      int i=0;
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it, ++i) {
	int mi=0;
	double sumH = 0.0, sumHenter = 0.0;
	vec gradi = zeros<vec>(n-1);
	vec gradHi = zeros<vec>(n-1);
	vec gradH0i = zeros<vec>(n-1);
	// vec grconstraint = zeros<vec>(n-1);
	for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	  sumH += H(*j);
	  gradHi += gradH.row(*j).t();
	  if (this->event(*j)==1) {
	    gradi += gradh.row(*j).t() / h(*j);
	    mi++;
	  }
	  if (this->delayed  && this->ind0[*j]) {
	    int jj = this->map0[*j];
	    sumHenter += H0(jj);
	    gradH0i += gradH0.row(jj).t();
	    grconstraint(i,span(0,n-2)) += this->kappa * (gradH0.row(jj) * H0(jj) * (H0(jj)<0));
	  }
	  grconstraint(i,span(0,n-2)) += this->kappa * ((gradh.row(*j) * h(*j) * (h(*j)<0)) + (gradH.row(*j) * H(*j) * (H(*j)<0)));
	}
	for (int k=0; k<n-1; ++k) {
	  if (recurrent) {
	    gr(i,k) += gradi(k) - (1.0/theta+mi)*theta*(gradHi(k)-gradH0i(k))/(1+theta*(sumH-sumHenter)) - grconstraint(k); // Rondeau et al (2005)
	  } else {
	    gr(i,k) += gradi(k) - theta*(1.0/theta+mi)*gradHi(k)/(1+theta*sumH) + 1.0/(1+theta*sumHenter)*gradH0i(k) - grconstraint(k); // Rondeau et al
	  }
	}
	if (recurrent) {
	  // axiom: D(-(1/exp(btheta)+mi)*log(1+exp(btheta)*H),btheta)
	  // axiom: D(log(Gamma(1/exp(btheta)+mi))-log(Gamma(1/exp(btheta)))+mi*btheta,btheta)
	  double H = sumH - sumHenter;
	  // gr(i,n-1) += ((1+theta*H)*log(1+H*theta)-H*mi*theta*theta - H*theta)/(H*theta*theta + theta) + (-R::digamma(1.0/theta+mi)+R::digamma(1/theta)+mi*theta)/theta;
	  gr(i,n-1) += ((1+theta*H)*log(1+H*theta)-H*mi*theta*theta - H*theta)/(H*theta*theta + theta);
	} else {
	  // axiom: D(-(1/exp(btheta)+mi)*log(1+exp(btheta)*H),btheta)
	  // axiom: D(1.0/exp(btheta)*log(1+exp(btheta)*H0),btheta)
	  gr(i,n-1) += ((1+sumH*theta)*log(1+sumH*theta)-sumH*mi*theta*theta-sumH*theta)/(sumH*theta*theta+theta);
	  gr(i,n-1) += (-(1+sumHenter*theta)*log(sumHenter*theta+1)+sumHenter*theta)/(sumHenter*theta*theta+theta); 
	}
	if (mi>0) {
		for (int k=1; k<=mi; ++k)
			// axiom: D(log(1+exp(btheta)*(mi-k)),btheta)
			gr(i,n-1) += theta*(mi-k)/(1.0+theta*(mi-k));
	}
      }
      gradli_constraint grli = {gr, grconstraint};
      return grli;  
    }
    bool feasible(vec beta) {
      vec coef(beta);
      coef.resize(this->n-1);
      return Base::feasible(coef);
    }
    OPTIM_FUNCTIONS;
    IndexMap clusters, cluster_delayed;
    bool recurrent;
  };


  /** @brief Wrapper for calling an objective_cluster method
   **/
  template<class T>
  double call_objective_cluster(double bi, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    return model->objective_cluster(bi);
  }    
  template<class T>
  double call_gradient_cluster(double bi, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    return model->gradient_cluster(bi);
  }    
  /** @brief Wrapper for calling a gradient_cluster method
   **/
  template<class T>
  double call_objective_cluster_bfgs(int n, double *bi, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    return model->objective_cluster(*bi);
  }    
  template<class T>
  void call_gradient_cluster_bfgs(int dim, double * bi, double* out, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    *out = model->gradient_cluster(*bi);
  }    
  /** 
      Extension of stpm2 and pstpm2 to include log-normal shared frailties with a cluster variable
  **/
  template<class Base>
  class NormalSharedFrailty : public Base {
  public:
    typedef std::vector<int> Index;
    typedef std::map<int,Index> IndexMap;
    typedef NormalSharedFrailty<Base> This;
    typedef std::map<int,double> Doubles;
    NormalSharedFrailty(SEXP sexp) : Base(sexp) {
      List list = as<List>(sexp);
      const BaseData fullTmp = {this->X, this->XD, this->X0, this->X1, this->bhazard,
				this->wt, this->wt0, this->event, this->time,
				this->map0, this->ind0, this->which0};
      full = fullTmp;
      IntegerVector cluster = as<IntegerVector>(list["cluster"]);
      full_Z = Z = as<vec>(list["Z"]);
      full_Z0 = Z0 = Z(this->ind0); // this could be size 0...
      gauss_x = as<vec>(list["gauss_x"]); // node values
      gauss_w = as<vec>(list["gauss_w"]); // probability weights
      adaptive = as<bool>(list["adaptive"]);
      recurrent = as<bool>(list["recurrent"]);
      first = true;
      // wragged array indexed by a map of vectors
      for (int id=0; id<cluster.size(); ++id) {
	clusters[cluster[id]].push_back(id);
      }
      this->nbeta = this->n - 1;
    }
    /** 
	Objective function for a log-normal shared frailty
	Assumes that weights == 1.
    **/
    double objective_nonadaptive(vec beta) {
      int n = beta.size();
      vec vbeta(beta); // nu=log(variance) is the last parameter in beta; exp(nu)=sigma^2
      vbeta.resize(n-1);
      // case: alpha for joint model?
      double sigma = exp(0.5*beta[n-1]); // standard deviation
      int K = gauss_x.size(); // number of nodes
      vec wstar = gauss_w/sqrt(datum::pi);
      double constraint = 0.0;
      double ll = 0.0;
      bool left_trunc_not_recurrent = (this->delayed && !this->recurrent && !this->interval);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	vec eta = this->X * vbeta;
	vec etaD = this->XD * vbeta;
	vec eta0(1,fill::zeros);
	vec eta1(this->X.n_rows,fill::zeros);
	if (this->delayed) {
	  eta0 = this->X0 * vbeta;
	  eta1 = this->X1 * vbeta;
	}
	double Lj = 0.0, L0j = 0.0;
	for (int k=0; k<K; ++k) {
	  li_constraint lik, H0ik;
	  double bi = sqrt(2.0)*sigma*this->gauss_x(k);
	  if (left_trunc_not_recurrent) {
	    this->delayed = false; 
	    lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	    this->delayed = true; 
	    H0ik = Base::li_left_truncated(eta0+Z0*bi);
	    L0j += exp(-sum(H0ik.li))*wstar(k);
	    constraint += H0ik.constraint;
	  } else {
	    lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	  }
	  Lj += exp(sum(lik.li))*wstar(k);
	  constraint += lik.constraint;
	}
	ll += log(Lj);
	if (left_trunc_not_recurrent)
	  ll -= log(L0j);
      }
      resetDesign();
      return -(ll - constraint);
    }
    double objective(vec beta) {
      return this->adaptive ? objective_adaptive(beta) : objective_nonadaptive(beta);
    }
    double objective_adaptive(vec beta) {
      int n = beta.size();
      this->objective_cluster_beta = beta; // for use by objective_cluster()
      vec vbeta(beta); // nu=log(variance) is the last parameter in beta; exp(nu)=sigma^2
      double sigma = exp(0.5*vbeta(this->n-1)); // standard deviation
      vbeta.resize(n-1);
      int K = gauss_x.size(); // number of nodes
      double constraint = 0.0;
      double ll = 0.0;
      bool left_trunc_not_recurrent = (this->delayed && !this->recurrent && !this->interval);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	// cluster-specific mode
	double mu;
	if (!first) {
	  // double Tol = this->reltol;
	  // int Maxit = 100;
	  // mu = R_zeroin2(-100.0,
	  // 		 100.0,
	  // 		 gradient_cluster(-100.0),
	  // 		 gradient_cluster(100.0),
	  // 		 &(call_gradient_cluster<This>),
	  // 		 (void *) this,
	  // 		 &Tol,
	  // 		 &Maxit);
	  mu = Brent_fmin(muhat[it->first]-1e-1,
	  		  muhat[it->first]+1e-1,
	  		  &(call_objective_cluster<This>),
	  		  (void *) this,
	  		  this->reltol);
	}
	if (first || std::abs(mu-muhat[it->first])>(1e-1-1e-5))
	// if (first || mu>100.0-1e-5 || mu<-100.0+1e-5)
	  mu = Brent_fmin(-100.0,100.0,&(call_objective_cluster<This>),(void *) this,
			  this->reltol);
	muhat[it->first] = mu;
	// cluster-specific variance
	double eps = 1e-6;
	double tau = sqrt(eps*2.0/(gradient_cluster(mu+eps)-gradient_cluster(mu-eps)));
	tauhat[it->first] = tau;
	vec eta = this->X * vbeta;
	vec etaD = this->XD * vbeta;
	vec eta0(1,fill::zeros);
	vec eta1(this->X.n_rows,fill::zeros);
	if (this->delayed) {
	  eta0 = this->X0 * vbeta;
	  eta1 = this->X1 * vbeta;
	}
	double Lj = 0.0, L0j = 0.0;
	for (int k=0; k<K; ++k) {
	  double g = 0.0, g0 = 0.0;
	  double bi = mu + sqrt(2.0)*tau*this->gauss_x(k);
	  li_constraint lik, l0ik;
	  if (left_trunc_not_recurrent) {
	    // first: reset such that delayed is true and calculate the likelihood for right censoring
	    this->delayed = false; 
	    lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	    // second: reset for delayed is true and calculate the likelihood for left truncation
	    this->delayed = true; 
	    l0ik = Base::li_left_truncated(eta0+Z0*bi);
	    constraint += l0ik.constraint;
	    g0 = exp(sum(l0ik.li)+R::dnorm(bi,0.0,sigma,1));
	    L0j += sqrt(2.0)*tau*g0*gauss_w(k)*exp(gauss_x(k)*gauss_x(k));
	  } else {
	    lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	  }
	  g = exp(sum(lik.li)+R::dnorm(bi,0.0,sigma,1));
	  Lj += sqrt(2.0)*tau*g*gauss_w(k)*exp(gauss_x(k)*gauss_x(k));
	  constraint += lik.constraint;
	}
	ll += log(Lj);
	if (left_trunc_not_recurrent && L0j>0.0)
	  ll -= log(L0j);
      }
      resetDesign();
      first = false;
      return -(ll - constraint);
    }
    vec getli(vec beta) { vec missing(1,fill::zeros); return missing; }
    mat getgradli(vec beta) { mat missing(1,1,fill::zeros); return missing; }
    // TODO: redefine li
    void resetDesign() {
      this->X = full.X;
      this->XD = full.XD;
      this->bhazard = full.bhazard;
      this->wt = full.wt;
      this->event = full.event;
      this->time = full.time;
      this->X1 = full.X1;
      this->X0 = full.X0;
      this->wt0 = full.wt0;
      this->which0 = full.which0;
      this->Z = full_Z;
      this->Z0 = full_Z0;
    }
    void clusterDesign(IndexMap::iterator it) {
      clusterid = it->first;
      uvec index = conv_to<uvec>::from(it->second);
      this->X = full.X.rows(index);
      this->XD = full.XD.rows(index);
      this->X1 = full.X1.rows(index);
      this->bhazard = full.bhazard(index);
      this->wt = full.wt(index);
      this->event = full.event(index);
      this->time = full.time(index);
      this->Z = full_Z(index);
      this->Z0 = vec(1,fill::zeros);
      if (this->delayed) {
	uvec index0 = Base::map0f(index);
	this->X0 = full.X0.rows(index0);
	this->wt0 = full.wt0(index0);
	this->Z0 = full_Z0(index0);
	this->which0 = Base::which0f(index);
      }
    }
    // log(g(bi))
    double objective_cluster(double bi) {
      vec vbeta = this->objective_cluster_beta; // nu=log(variance) is the last parameter in vbeta
      double sigma = exp(0.5*vbeta(this->n-1)); // standard deviation
      vbeta.resize(this->n - 1);
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec eta0(1,fill::zeros);
      vec eta1(this->X.n_rows,fill::zeros);
      if (this->delayed) {
	eta0 = this->X0 * vbeta;
	eta1 = this->X1 * vbeta;
      }
      li_constraint lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
      double ll = sum(lik.li) + R::dnorm(bi,0.0,sigma,1);
      return -ll;
    }
    double gradient_cluster(double bi) {
      vec vbeta = this->objective_cluster_beta; // nu=log(variance) is the last parameter in vbeta
      double sigma = exp(0.5*vbeta(this->n-1)); // standard deviation
      vbeta.resize(this->n - 1);
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec eta0(1,fill::zeros);
      vec eta1(this->X1.n_rows,fill::zeros);
      if (this->delayed) {
	eta0 = this->X0 * vbeta;
	eta1 = this->X1 * vbeta;
      }
      mat X = mat(Z); // mat(this->X.n_rows,1,fill::ones);
      mat XD = mat(this->XD.n_rows,1,fill::zeros);
      mat X0 = mat(Z0); // mat(this->X0.n_rows,1,fill::ones);
      mat X1 = mat(Z); // mat(this->X1.n_rows,1,fill::ones);
      gradli_constraint gradlik = Base::gradli(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi,X,XD,X0,X1);
      rowvec gradll = sum(gradlik.gradli,0) - bi/sigma/sigma;
      return -gradll(0);
    }
    // double fdhessian(boost::function<double(double x)> f, double x, double eps = 1.0e-6) {
    //   return (-f(x+2*eps)+16*f(x+eps)-30*f(x)+16*f(x-eps)-f(x-2*eps))/12.0/eps/eps;
    // }
    void calculate_modes_and_variances() {
      this->objective_cluster_beta = this->bfgs.coef;
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	modes[it->first] = Brent_fmin(-100.0,100.0,&(call_objective_cluster<This>),(void *) this,
				      this->reltol);
	// Abramowitz and Stegun 1972, p. 884
	double bi = modes[it->first];
	double eps = 1e-6;
	double f1 = objective_cluster(bi+2*eps);
	double f2 = objective_cluster(bi+eps);
	double f3 = objective_cluster(bi);
	double f4 = objective_cluster(bi-eps);
	double f5 = objective_cluster(bi-2*eps);
	double hessian = (-f1+16*f2-30*f3+16*f4-f5)/12.0/eps/eps;
	variances[it->first] = 1.0/hessian;
	variances[it->first] = eps*2.0/(gradient_cluster(bi+eps)-gradient_cluster(bi-eps));
      }
      resetDesign();
    }
    vec gradient(vec beta) {
      return this->adaptive ? gradient_adaptive(beta) : gradient_nonadaptive(beta);
    }
    vec gradient_adaptive(vec beta) {
      int n = beta.size();
      this->objective_cluster_beta = beta; // for use by objective_cluster()
      vec betastar = beta; // we use the last element as bi - is this a good idea?
      vec vbeta(beta); // nu=log(variance) is the last parameter in beta; exp(nu)=sigma^2
      double sigma = exp(0.5*vbeta(this->n-1)); // standard deviation
      vbeta.resize(n-1);
      int K = gauss_x.size(); // number of nodes
      // rowvec dconstraint(n,fill::zeros); 
      rowvec gradLi(n,fill::zeros);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	// cluster-specific modes and variances (calculated in the objective function)
	double mu = muhat[it->first];
	double tau = tauhat[it->first];
	mat Zmain(Z); 
	mat ZD(this->XD.n_rows,1,fill::zeros);
	mat Z0main(Z0); 
	mat Z1(Z); 
	mat Xstar = join_horiz(this->X, Zmain); 
	mat XDstar = join_horiz(this->XD, ZD); // assumes time-invariant random effects
	mat X0star = join_horiz(this->X0, Z0main); 
	mat X1star = join_horiz(this->X1, Z1);
	double Lj = 0.0;
	rowvec numerator(this->n,fill::zeros);
	for (int k=0; k<K; ++k) {
	  double bi = mu + sqrt(2.0)*tau*this->gauss_x(k);
	  betastar[betastar.size()-1] = bi;
	  vec etastar = Xstar * betastar;
	  vec etaDstar = XDstar * betastar;
	  vec eta0star = X0star * betastar;
	  vec eta1star = X1star * betastar;
	  li_constraint lik = Base::li(etastar,etaDstar,eta0star,eta1star);
	  double g = exp(sum(lik.li) + R::dnorm(bi,0.0,sigma,1));
	  gradli_constraint gradlik = Base::gradli(etastar, etaDstar, eta0star, eta1star,Xstar, XDstar, X0star, X1star);
	  Lj += sqrt(2.0)*tau*g*gauss_w(k)*exp(gauss_x(k)*gauss_x(k));
	  rowvec numeratorstar = sqrt(2.0)*tau*g*sum(gradlik.gradli,0)*gauss_w(k)*exp(gauss_x(k)*gauss_x(k));
	  numerator(span(0,n-2)) += numeratorstar(span(0,n-2));
	  numerator(n-1) += sqrt(2.0)*tau*gauss_w(k)*exp(gauss_x(k)*gauss_x(k))*g*(- 0.5 + bi*bi/2/sigma/sigma);
	  // dconstraint += sum(gradlik.constraint,0);
	}
	gradLi += numerator/Lj;
      }
      vec gr(n,fill::zeros);
      for (int i=0; i<n; i++) gr(i) = gradLi(i); // - dconstraint(i);
      resetDesign();
      return -gr;
    }
    vec gradient_nonadaptive(vec beta) {
	int n = beta.size();
	double sigma = exp(0.5*beta[n-1]);
	mat Zmain(Z);
	mat ZD(this->N,1,fill::zeros);
	mat Z0main(Z0);
	mat Z1(Z);
	mat Xstar = join_horiz(this->X, Zmain); 
	mat XDstar = join_horiz(this->XD, ZD); // assumes time-invariant random effects
	mat X0star = join_horiz(this->X0, Z0main); 
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
	// if (this->bfgs.trace>0) {
	// 	Rprintf("fdgradient="); Rprint(this->fdgradient(beta)); 
	// }
	vec gr(n,fill::zeros);
	for (int i=0; i<n; i++) gr(i) = grad(i); // - dconstraint(i);
	return -gr;  
    }
    bool feasible(vec beta) {
      vec coef(beta);
      coef.resize(beta.size()-1);
      return Base::feasible(coef);
    }
    SEXP return_modes() {
      calculate_modes_and_variances();
      return wrap(modes);
    }
    SEXP return_variances() {
      calculate_modes_and_variances();
      return wrap(variances);
    }
    OPTIM_FUNCTIONS;
    IndexMap clusters;
    vec gauss_x, gauss_w, full_Z, full_Z0, Z, Z0;
    BaseData full;
    Doubles modes, variances;
    vec objective_cluster_beta;
    Doubles muhat, tauhat;
    int clusterid;
    bool adaptive, first, recurrent;
  };

  template<class T>
  double call_objective_clusterND(int n, double * beta, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    vec bi(beta, n);
    return model->objective_cluster(bi);
  }    
  template<class T>
  void call_gradient_clusterND(int n, double * beta, double * out, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    vec bi(beta, n);
    vec grad = model->gradient_cluster(bi);
    for (int i=0; i<n; ++i) out[i] = grad[i];
  }    
  template<class Base>
  class NormalSharedFrailty2D : public Base {
  public:
    typedef std::vector<int> Index;
    typedef std::map<int,Index> IndexMap;
    typedef NormalSharedFrailty2D<Base> This;
    typedef std::map<int,vec> MapVec;
    typedef std::map<int,mat> MapMat;
    typedef std::map<int,cube> MapCube;
    NormalSharedFrailty2D(SEXP sexp) : Base(sexp) {
      List list = as<List>(sexp);
      const BaseData fullTmp = {this->X, this->XD, this->X0, this->X1, this->bhazard,
				this->wt, this->wt0, this->event, this->time,
				this->map0, this->ind0, this->which0};
      full = fullTmp;
      IntegerVector cluster = as<IntegerVector>(list["cluster"]);
      full_Z = Z = as<mat>(list["Z"]);
      full_Z0 = Z0 = Z.rows(this->ind0); // this could be size 0...
      gauss_x = as<vec>(list["gauss_x"]); // node values
      gauss_w = as<vec>(list["gauss_w"]); // probability weights
      adaptive = as<bool>(list["adaptive"]); // current = false
      recurrent = as<bool>(list["recurrent"]);
      first = true;
      // wragged array indexed by a map of vectors
      for (int id=0; id<cluster.size(); ++id) {
	clusters[cluster[id]].push_back(id);
      }
      redim = 2;
      reparm = 3;
      Sigma.set_size(redim,redim);
      invSigma.set_size(redim,redim);
      SqrtSigma.set_size(redim,redim);
      this->nbeta = this->n - reparm;
    }
    /** 
	Objective function for normal random effects
    **/
    double corrtransform(double x) { return 2.0/(1.0+exp(-x))-1.0; }
    double objective_nonadaptive(vec beta) {
      int n = beta.size();
      vec vbeta(beta);
      calc_SqrtSigma(beta);
      vbeta.resize(n-reparm);
      int K = gauss_x.size(); // number of nodes
      vec wstar = gauss_w/sqrt(datum::pi);
      vec d = sqrt(2.0)*gauss_x;
      double constraint = 0.0;
      double ll = 0.0;
      bool left_trunc_not_recurrent = (this->delayed && !this->recurrent && !this->interval);
      vec dk(redim), bi(redim);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	vec eta = this->X * vbeta;
	vec etaD = this->XD * vbeta;
	vec eta0(1,fill::zeros);
	vec eta1(this->X.n_rows,fill::zeros);
	if (this->delayed) {
	  eta0 = this->X0 * vbeta;
	  eta1 = this->X1 * vbeta;
	}
	double Lj = 0.0, L0j = 0.0;
	for (int k=0; k<K; ++k) {
	  for (int kk=0; kk<K; ++kk) {
	    li_constraint lik, H0ik;
	    dk(0) = d(k);
	    dk(1) = d(kk); 
	    bi = SqrtSigma * dk;
	    if (left_trunc_not_recurrent) {
	      this->delayed = false; 
	      lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	      this->delayed = true; 
	      H0ik = Base::li_left_truncated(eta0+Z0*bi);
	      L0j += exp(-sum(H0ik.li))*wstar(k)*wstar(kk);
	      constraint += H0ik.constraint;
	    } else {
	      lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	    }
	    Lj += exp(sum(lik.li))*wstar(k)*wstar(kk);
	    constraint += lik.constraint;
	  }
	}	  
	ll += log(Lj);
	if (left_trunc_not_recurrent)
	  ll -= log(L0j);
      }
      resetDesign();
      return -(ll - constraint);
    }
    double objective_cluster(vec bi) {
      vec vbeta = this->objective_cluster_beta; 
      vbeta.resize(this->n - reparm);
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec eta0(this->Z0.n_rows,fill::zeros);
      vec eta1(this->X.n_rows,fill::zeros);
      if (this->delayed) {
	eta0 = this->X0 * vbeta;
	eta1 = this->X1 * vbeta;
      }
      li_constraint lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
      vec zero(redim,fill::zeros);
      double ll = sum(lik.li) + dmvnrm_arma(bi,zero,this->Sigma,true);
      return -ll;
    }
    vec gradient_cluster(vec bi) {
      vec vbeta = this->objective_cluster_beta; 
      vbeta.resize(this->n - reparm);
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec eta0(1,fill::zeros);
      vec eta1(this->X1.n_rows,fill::zeros);
      if (this->delayed) {
	eta0 = this->X0 * vbeta;
	eta1 = this->X1 * vbeta;
      }
      mat X = Z; 
      mat XD = mat(this->XD.n_rows,redim,fill::zeros);
      mat X0 = Z0; 
      mat X1 = Z; 
      gradli_constraint gradlik = Base::gradli(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi,X,XD,X0,X1);
      vec gradll = (sum(gradlik.gradli,0)).t() - invSigma*bi;
      return -gradll;
    }
    mat calc_SqrtSigma(vec beta, bool calc_inverse = true) {
      bool flag;
      int n = beta.size();
      Sigma.resize(redim,redim);
      Sigma(0,0) = exp(beta[n-reparm]);       // variance
      double corr = corrtransform(beta[n-2]); // correlation
      Sigma(1,1) = exp(beta[n-1]);            // variance 
      Sigma(0,1) = Sigma(1,0) = corr*sqrt(Sigma(1,1)*Sigma(0,0)); // covariance
      if (calc_inverse)
	flag = inv(invSigma, Sigma);
      flag = chol(SqrtSigma,Sigma); // which version of RcppArmadillo??
      if (!flag) { Rprintf("Sigma:\n"); Rprint(Sigma); Rcpp::stop("Matrix sqrt invalid"); } 
      return SqrtSigma;
    }
    // gradient in SqrtSigma wrt beta
    cube gradSqrtSigma(vec beta, double eps = 1.0e-6) {
      cube out(redim,redim,reparm,fill::zeros);
      int offset = beta.size()-reparm;
      vec betax;
      mat val;
      for (int i=0; i<reparm; ++i) {
	betax = beta;
	betax(offset+i) += eps;
	val = calc_SqrtSigma(betax, false);
	betax(offset+i) -= 2.0*eps;
	out.slice(i) = (val - calc_SqrtSigma(betax, false))/2.0/eps;
      }
      return out;
    }
    // Sigma and SqrtSigma for the cluster-specific modes
    mat calc_SqrtSigma_adaptive(BFGS opt) {
      arma::mat tau, sqrttau;
      bool flag = inv(tau,
		      as<mat>(opt.calc_hessian(call_gradient_clusterND<This>, (void *) this)));
      flag = chol(sqrttau,tau); // which version of RcppArmadillo??
      if (!flag) { Rprintf("tau:\n"); Rprint(tau); Rcpp::stop("Matrix sqrt invalid"); }
      return sqrttau;
    }
    // gradient for SqrtSigma for the cluster-specific modes
    cube gradSqrtSigma_adaptive(BFGS opt, double eps = 1.0e-6) {
      cube out(redim,redim,reparm,fill::zeros);
      vec betax;
      vec old = this->objective_cluster_beta;
      mat val;
      for (int i=0; i<reparm; ++i) {
	betax = old;
	betax(i) += eps;
	this->objective_cluster_beta = betax;
	val = calc_SqrtSigma_adaptive(opt);
	betax(i) -= 2.0*eps;
	this->objective_cluster_beta = betax;
	out.slice(i) = (val - calc_SqrtSigma_adaptive(opt))/2.0/eps;
      }
      this->objective_cluster_beta = old;
      return out;
    }
    double objective_adaptive(vec beta) {
      int n = beta.size();
      this->objective_cluster_beta = beta; // for use by objective_cluster()
      vec vbeta(beta);
      calc_SqrtSigma(beta);
      vbeta.resize(n-reparm);
      int K = gauss_x.size(); // number of nodes
      vec wstar = gauss_w/sqrt(datum::pi);
      vec d = sqrt(2.0)*gauss_x;
      double constraint = 0.0;
      double ll = 0.0;
      bool left_trunc_not_recurrent = (this->delayed && !this->recurrent && !this->interval);
      vec dk(redim), bi(redim);
      vec zero(redim,fill::zeros);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	if (muhat.count(it->first)==0) { // initialise
	  muhat[it->first].set_size(redim);
	  for (int i=0; i<redim; ++i)
	    muhat[it->first](i) = 0.0;
	}
	// cluster-specific modes
	BFGS opt; 
	opt.optim(&(call_objective_clusterND<This>),
		  &(call_gradient_clusterND<This>),
		 as<NumericVector>(wrap(muhat[it->first])),
		 (void *) this);
	muhat[it->first] = as<vec>(wrap(opt.coef));
	// variance for the cluster-specific modes
	arma::mat tau, sqrttau;
	bool flag = inv(tau,
			as<mat>(opt.calc_hessian(call_gradient_clusterND<This>, (void *) this)));
	double dettau = det(tau); // logdet would be more precise
	dettauhat[it->first] = dettau;
	flag = chol(sqrttau,tau); // which version of RcppArmadillo?
	if (!flag) { Rprintf("tau:\n"); Rprint(tau); Rcpp::stop("Matrix sqrt invalid."); }
	sqrttauhat[it->first] = sqrttau;
	if (this->optimiser=="BFGS")
	  gradsqrttauhat[it->first] = gradSqrtSigma_adaptive(opt); 
	vec eta = this->X * vbeta;
	vec etaD = this->XD * vbeta;
	vec eta0(1,fill::zeros);
	vec eta1(this->X.n_rows,fill::zeros);
	if (this->delayed) {
	  eta0 = this->X0 * vbeta;
	  eta1 = this->X1 * vbeta;
	}
	double Lj = 0.0, L0j = 0.0;
	for (int k=0; k<K; ++k) {
	  for (int kk=0; kk<K; ++kk) {
	    dk(0) = d(k);
	    dk(1) = d(kk); 
	    bi = muhat[it->first] + sqrttau * dk;
	    double g = 0.0, g0 = 0.0;
	    li_constraint lik, l0ik;
	    if (left_trunc_not_recurrent) {
	      this->delayed = false; 
	      lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	      this->delayed = true; 
	      l0ik = Base::li_left_truncated(eta0+Z0*bi);
	      g0 = exp(sum(l0ik.li)+dmvnrm_arma(bi,zero,Sigma,true));
	      L0j += wstar(k)*wstar(kk)*g0*exp(dot(d,d)/2);
	    } else {
	      lik = Base::li(eta+Z*bi,etaD,eta0+Z0*bi,eta1+Z*bi);
	    }
	    g = exp(sum(lik.li)+dmvnrm_arma(bi,zero,Sigma,true));
	    Lj += wstar(k)*wstar(kk)*g*exp(dot(d,d)/2);
	    constraint += lik.constraint;
	  }
	}
	ll += log(2.0*datum::pi)*redim/2.0+log(dettau)/2.0+log(Lj);
	if (left_trunc_not_recurrent)
	  ll -= log(L0j);
      }
      resetDesign();
      first = false;
      return -(ll - constraint);
    }
    vec gradient(vec beta) {
      return this->adaptive ? gradient_adaptive(beta) : gradient_nonadaptive(beta);
    }
    vec gradient_nonadaptive(vec beta) {
      int n = beta.size();
      vec vbeta(beta);
      vbeta.resize(n-reparm);
      calc_SqrtSigma(beta);
      cube gradSqrt = gradSqrtSigma(beta);
      if (this->bfgs.trace>1)
	Rprint(gradSqrt);
      int K = gauss_x.size(); // number of nodes
      vec wstar = gauss_w/sqrt(datum::pi);
      vec d = sqrt(2.0)*gauss_x;
      vec gradl(n,fill::zeros);
      vec dk(redim), bi(redim);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	uvec index = conv_to<uvec>::from(it->second);
	mat Zmain(Z);
	mat ZD(this->XD.n_rows,redim,fill::zeros);
	mat Z0main(Z0);
	mat Z1(Z);
	mat Xstar = join_horiz(this->X, Zmain); 
	mat XDstar = join_horiz(this->XD, ZD); // assumes time-invariant random effects
	mat X0star = join_horiz(this->X0, Z0main); 
	mat X1star = join_horiz(this->X1, Z1);
	double L = 0.0;
	vec gradL(n,fill::zeros);
	vec betastar = beta; betastar.resize(beta.size() - reparm + redim);
	// for K nodes calculation 
	for (int k=0; k<K; ++k) {
	  for (int kk=0; kk<K; ++kk) {
	    dk(0) = d(k);
	    dk(1) = d(kk); 
	    bi = SqrtSigma * dk;
	    cube dbic = gradSqrt.each_slice() % dk;
	    mat dbi = reshape( mat(dbic.memptr(), dbic.n_elem, 1, false), redim, reparm);
	    betastar(span(betastar.size()-redim,betastar.size()-1)) = bi;
	    vec etastar = Xstar * betastar;
	    vec etaDstar = XDstar * betastar;
	    vec eta0star = X0star * betastar;
	    vec eta1star = X1star * betastar;
	    li_constraint lik = Base::li(etastar, etaDstar, eta0star, eta1star);
	    gradli_constraint gradlik = Base::gradli(etastar, etaDstar, eta0star, eta1star,Xstar, XDstar, X0star, X1star);
	    // adjust the last columns of the gradient to account for the variance components:
	    // chain rule: d lik/d bi * d bi/d nu
	    int restart = gradlik.gradli.n_cols-redim, restop = gradlik.gradli.n_cols-1;
	    gradlik.gradli.resize(gradlik.gradli.n_rows, gradlik.gradli.n_cols-redim+reparm);
	    gradlik.gradli.cols(gradlik.gradli.n_cols-reparm,gradlik.gradli.n_cols-1) = gradlik.gradli.cols(restart, restop) * dbi;
	    double Lk = exp(sum(lik.li));
	    L += Lk*wstar(k)*wstar(kk);
	    gradL += (Lk*sum(gradlik.gradli,0)*wstar(k)*wstar(kk)).t();
	  }
	}
	gradl += gradL/L;
      }
      resetDesign();
      return -gradl;  
    }
    vec gradient_adaptive(vec beta) {
      int n = beta.size();
      objective_adaptive(beta); // bug check - remove
      this->objective_cluster_beta = beta; // for use by objective_cluster()
      vec vbeta(beta);
      vbeta.resize(n-reparm);
      calc_SqrtSigma(beta); // unnecessary?
      int K = gauss_x.size(); // number of nodes
      vec wstar = gauss_w/sqrt(datum::pi);
      vec d = sqrt(2.0)*gauss_x;
      vec gradl(n,fill::zeros);
      vec dk(redim), bi(redim);
      vec zero(redim,fill::zeros);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	uvec index = conv_to<uvec>::from(it->second);
	vec mu = muhat[it->first];
	mat sqrttau = sqrttauhat[it->first];
	cube gradsqrttau = gradsqrttauhat[it->first];
	mat Zmain(Z);
	mat ZD(this->XD.n_rows,redim,fill::zeros);
	mat Z0main(Z0);
	mat Z1(Z);
	mat Xstar = join_horiz(this->X, Zmain); 
	mat XDstar = join_horiz(this->XD, ZD); // assumes time-invariant random effects
	mat X0star = join_horiz(this->X0, Z0main); 
	mat X1star = join_horiz(this->X1, Z1);
	double L = 0.0;
	vec gradL(n,fill::zeros);
	vec betastar = beta; betastar.resize(beta.size() - reparm + redim);
	// for K nodes calculation 
	for (int k=0; k<K; ++k) {
	  for (int kk=0; kk<K; ++kk) {
	    dk(0) = d(k);
	    dk(1) = d(kk); 
	    bi = mu + sqrttau * dk;
	    cube dbic = gradsqrttau.each_slice() % dk;
	    mat dbi = reshape( mat(dbic.memptr(), dbic.n_elem, 1, false), redim, reparm);
	    betastar(span(betastar.size()-redim,betastar.size()-1)) = bi;
	    vec etastar = Xstar * betastar;
	    vec etaDstar = XDstar * betastar;
	    vec eta0star = X0star * betastar;
	    vec eta1star = X1star * betastar;
	    li_constraint lik = Base::li(etastar, etaDstar, eta0star, eta1star);
	    gradli_constraint gradlik = Base::gradli(etastar, etaDstar, eta0star, eta1star,Xstar, XDstar, X0star, X1star);
	    // adjust the last columns of the gradient to account for the variance components:
	    // chain rule: d lik/d bi * d bi/d nu
	    int restart = gradlik.gradli.n_cols-redim, restop = gradlik.gradli.n_cols-1;
	    gradlik.gradli.resize(gradlik.gradli.n_rows, gradlik.gradli.n_cols-redim+reparm);
	    gradlik.gradli.cols(gradlik.gradli.n_cols-reparm,gradlik.gradli.n_cols-1) = gradlik.gradli.cols(restart, restop) * dbi;
	    double Lk = exp(sum(lik.li)+dmvnrm_arma(bi,zero,Sigma,true));
	    L += Lk*wstar(k)*wstar(kk)*exp(dot(d,d)/2);
	    gradL += (Lk*sum(gradlik.gradli,0)*wstar(k)*wstar(kk)*exp(dot(d,d)/2)).t();
	  }
	}
	// missing: grad(log(det(tau))/2,beta)
	gradl += gradL/L;
      }
      resetDesign();
      return -gradl;  
    }
    double objective(vec beta) {
      return adaptive ? objective_adaptive(beta) : objective_nonadaptive(beta);
    }
    void resetDesign() {
      this->X = full.X;
      this->XD = full.XD;
      this->bhazard = full.bhazard;
      this->wt = full.wt;
      this->event = full.event;
      this->time = full.time;
      this->X1 = full.X1;
      this->X0 = full.X0;
      this->wt0 = full.wt0;
      this->which0 = full.which0;
      this->Z = full_Z;
      this->Z0 = full_Z0;
    }
    void clusterDesign(IndexMap::iterator it) {
      clusterid = it->first;
      uvec index = conv_to<uvec>::from(it->second);
      this->X = full.X.rows(index);
      this->XD = full.XD.rows(index);
      this->X1 = full.X1.rows(index);
      this->bhazard = full.bhazard(index);
      this->wt = full.wt(index);
      this->event = full.event(index);
      this->time = full.time(index);
      this->Z = full_Z.rows(index);
      this->Z0 = mat(1,redim,fill::zeros);
      if (this->delayed) {
	uvec index0 = Base::map0f(index);
	this->X0 = full.X0.rows(index0);
	this->wt0 = full.wt0(index0);
	this->Z0 = full_Z0.rows(index0);
	this->which0 = Base::which0f(index);
      }
    }
    bool feasible(vec beta) {
      vec coef(beta);
      coef.resize(beta.size()-reparm);
      return Base::feasible(coef);
    }
    void calculate_modes_and_variances() {
      vec beta = as<vec>(this->bfgs.coef);
      this->objective_cluster_beta = beta;
      calc_SqrtSigma(beta);
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	clusterDesign(it);
	if (muhat.count(it->first)==0) { // initialise
	  muhat[it->first].set_size(redim);
	  for (int i=0; i<redim; ++i)
	    muhat[it->first](i) = 0.0;
	}
	// cluster-specific modes
	BFGS opt; 
	opt.optim(&(call_objective_clusterND<This>),
		  &(call_gradient_clusterND<This>),
		  as<NumericVector>(wrap(muhat[it->first])),
		  (void *) this);
	muhat[it->first] = as<vec>(wrap(opt.coef));
	// cluster-specific variance
	arma::mat tau, sqrttau;
	bool flag = inv(tau,
			as<mat>(opt.calc_hessian(call_gradient_clusterND<This>, (void *) this)));
	double dettau = det(tau);
	dettauhat[it->first] = dettau;
	flag = chol(sqrttau,tau); // which version of RcppArmadillo?
	if (!flag) { Rprintf("tau:\n"); Rprint(tau); Rcpp::stop("Matrix sqrt invalid."); }
	sqrttauhat[it->first] = sqrttau;
	gradsqrttauhat[it->first] = gradSqrtSigma_adaptive(opt); 
      }
      resetDesign();
    }
    SEXP return_modes() {
      calculate_modes_and_variances();
      return wrap(muhat);
    }
    SEXP return_variances() { 
      calculate_modes_and_variances();
      return wrap(sqrttauhat); // actually, these are not variances!
    }
    OPTIM_FUNCTIONS;
    IndexMap clusters;
    vec gauss_x, gauss_w;
    mat full_Z, full_Z0, Z, Z0;
    BaseData full;
    MapVec modes;
    MapMat variances;
    vec objective_cluster_beta;
    MapVec muhat;
    MapMat sqrttauhat;
    MapCube gradsqrttauhat;
    std::map<int,double> dettauhat;
    int clusterid, redim, reparm;
    bool adaptive, first, recurrent;
    mat Sigma, SqrtSigma, invSigma;
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
      // model.bfgs.hessian = model.bfgs.calc_hessian(&optimgradient<T>, (void *) &model);
      // model.bfgs.hessian = model.calcHessian(model.bfgs.coef);
      model.post_process();
      return List::create(_("fail")=model.bfgs.fail, 
			  _("coef")=wrap(model.bfgs.coef),
			  _("hessian")=wrap(model.bfgs.hessian),
			  _("kappa")=wrap(model.kappa));
    }
    else if (return_type == "objective")
      return wrap(model.objective(beta));
    else if (return_type == "gradient") {
      model.objective(beta); // only needed for updating NormalSharedFrailty
      return wrap(model.gradient(beta));
    }
    else if (return_type == "feasible")
      return wrap(model.feasible(beta));
    else if (return_type == "modes")
      return model.return_modes();
    else if (return_type == "variances")
      return model.return_variances();
    else if (return_type == "li")
      return wrap(model.getli(beta));
    else if (return_type == "gradli")
      return wrap(model.getgradli(beta));
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
    else if (return_type == "gradient") {
      model.objective(beta);
      return wrap(model.gradient(beta));
    }
    else if (return_type == "gradient0")
      return wrap(model.gradient0(beta));
    else if (return_type == "constraint")
      return wrap(model.feasible(beta));
    else if (return_type == "feasible")
      return wrap(model.feasible(beta));
    else if (return_type == "li")
      return wrap(model.getli(beta));
    else if (return_type == "gradli")
      return wrap(model.getgradli(beta));
    else if (return_type == "modes")
      return model.return_modes();
    else if (return_type == "variances")
      return model.return_variances();
    else {
      REprintf("Unknown return_type.\n");
      return wrap(-1);
    }
  }
  RcppExport SEXP model_output(SEXP args) {
    List list = as<List>(args);
    std::string type = as<std::string>(list["type"]);
    if (type=="stpm2") {
      return stpm2_model_output_<Stpm2>(args);
    }
    else if (type=="pstpm2") {
      return pstpm2_model_output_<Pstpm2<Stpm2,SmoothLogH> >(args);
    }
    else if (type=="stpm2_gamma_frailty") {
      return stpm2_model_output_<GammaSharedFrailty<Stpm2> >(args);
    }
    else if (type=="pstpm2_gamma_frailty") {
      return pstpm2_model_output_<Pstpm2<GammaSharedFrailty<Stpm2>,SmoothLogH> >(args);
    }
    else if (type=="stpm2_normal_frailty") {
      return stpm2_model_output_<NormalSharedFrailty<Stpm2> >(args);
    }
    else if (type=="stpm2_normal_frailty_2d") {
      return stpm2_model_output_<NormalSharedFrailty2D<Stpm2> >(args);
    }
    else if (type=="pstpm2_normal_frailty") {
      // order is important here
      return pstpm2_model_output_<Pstpm2<NormalSharedFrailty<Stpm2>,SmoothLogH> >(args);
    }
    else if (type=="pstpm2_normal_frailty_2d") {
      // order is important here
      return pstpm2_model_output_<Pstpm2<NormalSharedFrailty2D<Stpm2>,SmoothLogH> >(args);
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
} // rstpm2 namespace
