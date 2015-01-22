#include <RcppArmadillo.h>
#include <vector>
#include <map>
#include "c_optim.h"

namespace rstpm2 {

  using namespace Rcpp;
  using namespace arma;

  // typedefs

  typedef bool constraintfn(int, double *, void *);

  // Hadamard element-wise multiplication for the _columns_ of a matrix with a vector

  mat rmult(mat m, vec v) {
    mat out(m);
    out.each_col() %= v;
    return out;
  }
  mat lmult(vec v, mat m) {
    return rmult(m,v);
  }

  // element-wise multiplication and division for NumericVector and arma::vec

  NumericVector ew_mult(NumericVector nv, vec v) {
    return as<NumericVector>(wrap(as<vec>(wrap(nv)) % v));
  }

  NumericVector ew_div(NumericVector nv, vec v) {
    return as<NumericVector>(wrap(as<vec>(wrap(nv)) / v));
  }

  void Rprint(NumericMatrix m) {
    for (int i=0; i<m.nrow(); ++i) {
      for (int j=0; j<m.ncol(); ++j) 
	Rprintf("%f ", m(i,j));
      Rprintf("\n");
    }
  }

  void Rprint(mat m) {
    for (size_t i=0; i<m.n_rows; ++i) {
      for (size_t j=0; j<m.n_cols; ++j) 
	Rprintf("%f ", m(i,j));
      Rprintf("\n");
    }
  }

  vec pnorm01(vec x) {
    vec out(x.size());
    for (size_t i=0; i<x.size(); ++i)
      out(i) = R::pnorm(x(i),0.0,1.0,1,0);
    return out;
    //return as<vec>(wrap(pnorm(as<NumericVector>(wrap<vec>(x)),0.0,1.0)));
  }

  vec qnorm01(vec x) {
    vec out(x.size());
    for (size_t i=0; i<x.size(); ++i)
      out(i) = R::qnorm(x(i),0.0,1.0,1,0);
    return out;
  }

  vec dnorm01(vec x) {
    vec out(x.size());
    for (size_t i=0; i<x.size(); ++i)
      out(i) = R::dnorm(x(i),0.0,1.0,0);
    return out;
  }

  vec logit(vec p) {
    return log(p/(1-p));
  }

  vec expit(vec x) {
    return 1/(1+exp(-x));
  }

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

  // Q: do we need parscale, trace, reltol and criterion in these structs?

  class stpm2 {
  public:
    stpm2(SEXP sexp) {
      List list = as<List>(sexp);
      init = as<NumericVector>(list["init"]);
      X = as<mat>(list["X"]); 
      XD = as<mat>(list["XD"]); 
      bhazard = as<vec>(list["bhazard"]);
      wt = as<vec>(list["wt"]);
      event = as<vec>(list["event"]);
      delayed = as<int>(list["delayed"]);
      X0.set_size(1,1); X0.fill(0.0);
      wt0.set_size(1); wt0.fill(0.0);
      if (delayed == 1) {
	X0 = as<mat>(list["X0"]); // optional
	XD0 = as<mat>(list["XD0"]); // optional
	wt0 = as<vec>(list["wt0"]); // optional
      }
      parscale = as<vec>(list["parscale"]);
      reltol = as<double>(list["reltol"]);
      kappa = as<double>(list["kappa"]);
      trace = as<int>(list["trace"]);
    }
    // defaults are for the PH model
    virtual vec link(vec S) { return log(-log(S)); }
    virtual vec ilink(vec x) { return exp(-exp(x)); }
    virtual vec h(vec eta, vec etaD) { return etaD % exp(eta); }
    virtual vec H(vec eta, vec etaD) { return exp(eta); }
    virtual mat gradH(vec eta, vec etaD, mat X, mat XD) { return rmult(X,exp(eta)); }
    virtual mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(XD, exp(eta)) + rmult(X, etaD % exp(eta));
    }
    NumericVector init;
    mat X, XD, X0, XD0; 
    vec bhazard,wt,wt0,event,parscale;
    double reltol, kappa;
    int delayed, trace;
  };

  class stpm2PO : public stpm2 {
  public:
    stpm2PO(SEXP sexp) : stpm2(sexp) { }
    vec link(vec S) { return -logit(S); }
    vec ilink(vec x) { return expit(-x); }
    vec h(vec eta, vec etaD) { return etaD % exp(eta) % expit(-eta); }
    vec H(vec eta, vec etaD) { return log(1+exp(eta)); }
    mat gradH(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(X,exp(eta) % expit(-eta));
    }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(X, etaD % exp(eta) % expit(-eta)) - 
	rmult(X, exp(2*eta) % etaD % expit(-eta) % expit(-eta)) +
	rmult(XD, exp(eta) % expit(-eta));
    }
  };

  class stpm2Probit : public stpm2 {
  public:
    stpm2Probit(SEXP sexp) : stpm2(sexp) { }
    vec link(vec S) { return -qnorm01(S); }
    vec ilink(vec eta) { return pnorm01(-eta); }
    vec H(vec eta, vec etaD) { return -log(pnorm01(-eta)); }
    vec h(vec eta, vec etaD) { return etaD % dnorm01(eta) / pnorm01(-eta); }
    mat gradH(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(X, dnorm01(eta) / pnorm01(-eta));
    }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(X, -eta % dnorm01(eta) % etaD / pnorm01(-eta)) +
	rmult(X, dnorm01(eta) % dnorm01(eta) / pnorm01(-eta) / pnorm01(-eta) % etaD) +
	rmult(XD,dnorm01(eta) / pnorm01(-eta));
    }
  };
  
  class stpm2AH : public stpm2 {
  public:
    stpm2AH(SEXP sexp) : stpm2(sexp) { }
    vec link(vec S) { return -log(S); }
    vec ilink(vec x) { return exp(-x); }
    vec h(vec eta, vec etaD) { return etaD; }
    vec H(vec eta, vec etaD) { return eta; }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { return XD; }
    mat gradH(vec eta, vec etaD, mat X, mat XD) { return X;  }
  };

  struct SmoothLogH {
    int first_para, last_para;
    mat S;
  };

  struct SmoothHaz {
    mat X0, X1, X2, X3;
    vec w;
    double lambda;
  };

  template<class Smooth>
  std::vector<Smooth> read_smoothers(List lsmooth) {}

  template<>
  std::vector<SmoothLogH> read_smoothers(List lsmooth) {
    std::vector<SmoothLogH> smooth;
    for(int i=0; i<lsmooth.size(); ++i) {
      List lsmoothi = as<List>(lsmooth[i]);
      List lsmoothiS = as<List>(lsmoothi["S"]);
      SmoothLogH smoothi =  {
	as<int>(lsmoothi["first.para"]) - 1, 
	as<int>(lsmoothi["last.para"]) - 1, 
	as<mat>(lsmoothiS[0])
      };
      smooth.push_back(smoothi);
    }
    return smooth;
  }

  template<>
  std::vector<SmoothHaz> read_smoothers(List lsmooth) {
    std::vector<SmoothHaz> smooth;
    for(int i=0; i<lsmooth.size(); ++i) {
      List lsmoothi = as<List>(lsmooth[i]);
      SmoothHaz smoothi =  {
	as<mat>(lsmoothi["X0"]), 
	as<mat>(lsmoothi["X1"]), 
	as<mat>(lsmoothi["X2"]), 
	as<mat>(lsmoothi["X3"]), 
	as<vec>(lsmoothi["w"]), 
	as<double>(lsmoothi["lambda"])
      };
      smooth.push_back(smoothi);
    }
    return smooth;
  }

  template<class Smooth, class Stpm2 = stpm2>
  class pstpm2 : public Stpm2 {
  public:
    pstpm2(SEXP sexp) : Stpm2(sexp) {
      List list = as<List>(sexp);
      List lsmooth = as<List>(list["smooth"]);
      sp = as<NumericVector>(list["sp"]);
      reltol_search = as<double>(list["reltol_search"]);
      // reltol = as<double>(list["reltol"]);
      alpha = as<double>(list["alpha"]);
      criterion = as<int>(list["criterion"]);
      smooth = read_smoothers<Smooth>(lsmooth);
    }
    NumericVector sp;
    double alpha, reltol_search; 
    int criterion;
    std::vector<Smooth> smooth;
  };

  /** 
      Extension of stpm2 and pstpm2 to include shared frailties with a cluster variable
  **/
  template<class Stpm2>
  class stpm2Shared : public Stpm2 {
  public:
    stpm2Shared(SEXP sexp) : Stpm2(sexp) {
      List list = as<List>(sexp);
      IntegerVector cluster = as<IntegerVector>(list["cluster"]);
      // wragged array indexed by a map of vectors
      for (int i=0; i<cluster.size(); ++i) {
	clusters[cluster[i]].push_back(i);
      }
    }
    std::map<int,std::vector<int> > clusters;
  };


  /**
     Extension to the BFGS class

     optim() and calc_hessian() assume that data->parscale is an array/vector/NumericVector
     optimWithConstraint() further assumes that data->kappa is a double with an initial value
   */
  template<class Data>
  class BFGS2 : public BFGS {
  public:
    void optim(optimfn fn, optimgr gr, NumericVector init, void * ex,
	       bool apply_parscale = true, double eps = 1.0e-8) {
      Data * data = (Data *) ex;
      n = init.size();
      if (apply_parscale) for (int i = 0; i<n; ++i) init[i] /= data->parscale[i];
      BFGS::optim(fn,gr,init,ex,eps);
      if (apply_parscale) for (int i = 0; i<n; ++i) coef[i] *= data->parscale[i];
      hessian = calc_hessian(gr, ex, eps);
    }
    void optimWithConstraint(optimfn fn, optimgr gr, NumericVector init, void * ex, constraintfn constraint,
			     bool apply_parscale = true, double eps = 1.0e-8) {
      Data * data = (Data *) ex;
      n = init.size();
      if (apply_parscale) for (int i = 0; i<n; ++i) init[i] /= data->parscale[i];
      bool satisfied;
      do {
	BFGS::optim(fn,gr,init,ex,eps);
	satisfied = constraint(n,&coef[0],ex);
	if (!satisfied) data->kappa *= 2.0;   
      } while ((!satisfied)&& data->kappa < 4*1.0e18);      
      if (apply_parscale) for (int i = 0; i<n; ++i) coef[i] *= data->parscale[i];
      hessian = calc_hessian(gr, ex, eps);
   // Rprintf("kappa=%f\n",data->kappa);
    }
    NumericMatrix calc_hessian(optimgr gr, void * ex, double eps = 1.0e-8) {
      Data * data = (Data *) ex;
      vec parscale(n);
      for (int i=0; i<n; ++i) {
	parscale[i] = data->parscale[i];
	data->parscale[i] = 1.0;
      }
      NumericMatrix hessian = BFGS::calc_hessian(gr,ex,eps);
      for (int i=0; i<n; ++i) data->parscale[i] = parscale[i];
      return hessian;
    }
  };

  template<class Data>
  double fminfn(int n, double * beta, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec h = data->h(eta,etaD) + data->bhazard;
    vec H = data->H(eta,etaD) + data->bhazard;
    double constraint = data->kappa/2.0 * sum(h % h % (h<0)); // sum(h^2 | h<0)
    vec eps = h*0.0 + 1.0e-16; 
    h = max(h,eps);
    double ll = sum(data->wt % data->event % log(h)) - sum(data->wt % H) - constraint;
    if (data->delayed == 1) {
      vec eta0 = data->X0 * vbeta;
      vec etaD0 = data->XD0 * vbeta;
      mat H0 = data->H(eta0,etaD0);
      ll += sum(data->wt0 % H0);
    }
    return -ll;  
  }

  /** 
      Objective function for a Gamma shared frailty
      Assumes that weights == 1.
      The implementation is incomplete, requiring either the use of optimisation without
      gradients or the calculation of the gradients.
  **/
  template<class Data>
  double fminfnShared(int n, double * beta, void *ex) {
    typedef std::vector<int> Index;
    typedef std::map<int,Index> IndexMap;
    Data * data = (Data *) ex;
    vec vbeta(beta,n-1); // theta is the last parameter in beta
    double theta = beta[n-1];
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec h = data->h(eta,etaD) + data->bhazard;
    vec H = data->H(eta,etaD) + data->bhazard;
    double constraint = data->kappa/2.0 * sum(h % h % (h<0)); // sum(h^2 | h<0)
    vec eps = h*0.0 + 1.0e-16; 
    h = max(h,eps);
    mat H0;
    //double ll = sum(data->wt % data->event % log(h)) - sum(data->wt % H) - constraint;
    if (data->delayed == 1) {
      vec eta0 = data->X0 * vbeta;
      vec etaD0 = data->XD0 * vbeta;
      H0 = data->H(eta0,etaD0);
      //ll += sum(data->wt0 % H0);
    }
    double ll;
    for (IndexMap::iterator it=data->clusters.begin(); it!=data->clusters.end(); ++it) {
      int mi=0;
      double sumH = 0.0, sumHenter = 0.0;
      for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	if (data->event(*j)==1) {
	  mi++;
	  ll += log(h(*j));
	}
	sumH += H(*j);
	if (data->delayed == 1)
	  sumHenter += H0(*j);
      }
      ll -= (1.0/theta+mi)*log(1.0+theta*sumH);
      if (data->delayed == 1)
	ll += 1.0/theta*log(1.0+theta*sumHenter);
      if (mi>0) {
	int k=1;
	for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	  if (data->event(*j)==1) {
	    ll += log(1.0+theta*(mi-k));
	    k++;
	  }
	}
      }
    }
    ll -= constraint;
    return -ll;  
  }

  /** 
      Utility function for calculating a gamma frailty log-likelihood within R code 
   **/
  RcppExport SEXP llgammafrailty(SEXP _theta, SEXP _h, SEXP _H, SEXP _d,
				    SEXP _cluster) {
    double theta = as<double>(_theta);
    NumericVector h = as<NumericVector>(_h);
    NumericVector H = as<NumericVector>(_H);
    IntegerVector d = as<IntegerVector>(_d);
    IntegerVector cluster = as<IntegerVector>(_cluster);
    double ll = 0.0;
    // wragged array indexed by a map of vectors
    typedef std::vector<int> Index;
    typedef std::map<int,Index> IndexMap;
    IndexMap clusters;
    for (int i=0; i<cluster.size(); ++i) {
      clusters[cluster[i]].push_back(i);
    }
    for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
      int mi=0;
      double sumH = 0.0;
      for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	if (d[*j]==1) {
	  mi++;
	  ll += log(h[*j]);
	}
	sumH += H[*j];
	// sumHenter += Henter[*j];
      }
      ll -= (1.0/theta+mi)*log(1.0+theta*sumH);
      // ll += 1.0/theta*log(1.0+theta*sumHenter);
      if (mi>0) {
	int k=1;
	for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	  if (d[*j]==1) {
	    ll += log(1.0+theta*(mi-k));
	    k++;
	  }
	}
      }
    }
    return wrap(ll);
  }


  template<class Data>
  bool fminfn_constraint(int n, double * beta, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec h = data->h(eta, etaD) + data->bhazard;
    return all(h>0);
  }

  template<class Data>
  void fminfn_nlm(int n, double * beta, double * f, void *ex) {
    *f = fminfn<Data>(n, beta, ex);
  };


  // log H penalty - using penalty matrices from mgcv
  template<class Stpm2>
  double pfminfn_SmoothLogH(int n, double * beta, void *ex) {
    typedef SmoothLogH Smooth;
    typedef pstpm2<Smooth,Stpm2> Data;
    double ll = -fminfn<Data>(n,beta,ex);
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    for (size_t i=0; i < data->smooth.size(); ++i) {
      Smooth smoothi = data->smooth[i];
      ll -= (data->sp)[i]/2 * 
	dot(vbeta.subvec(smoothi.first_para,smoothi.last_para),
	    smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para));
    }
    return -ll;  
  }

  // hazard penalty
  template<class Stpm2>
  double pfminfn_SmoothHaz(int n, double * beta, void *ex) {
    typedef SmoothHaz Smooth;
    typedef pstpm2<Smooth,Stpm2> Data;
    double ll = -fminfn<Data>(n,beta,ex);
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    for (size_t i=0; i < data->smooth.size(); ++i) {
      Smooth obj = data->smooth[i];
      vec s0 = obj.X0 * vbeta;
      vec s1 = obj.X1 * vbeta;
      vec s2 = obj.X2 * vbeta;
      vec s3 = obj.X3 * vbeta;
      vec h2 = exp(s0) % (s3 + 3*s1%s2 + s1%s1%s1);
      ll -= (data->sp)[i]/2 * obj.lambda * sum(obj.w % h2 % h2);
    }
    return -ll;  
  };

  template<class Data>
  void grfn(int n, double * beta, double * gr, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec h = data->h(eta,etaD) + data->bhazard;
    // vec H = exp(eta);
    mat gradH = data->gradH(eta,etaD,data->X,data->XD);
    mat gradh = data->gradh(eta,etaD,data->X,data->XD);
    mat Xgrad = -gradH + rmult(gradh, data->event / h);
    mat Xconstraint = - data->kappa*rmult(gradh,h);
    vec eps = h*0.0 + 1.0e-16; // hack
    //Xgrad.rows(h<=eps) = Xconstraint.rows(h<=eps);
    Xgrad = rmult(Xgrad,data->wt);
    rowvec vgr = sum(Xgrad,0);
    if (data->delayed == 1) {
      vec eta0 = data->X0 * vbeta;
      vec etaD0 = data->XD0 * vbeta;
      mat gradH0 = data->gradH(eta0, etaD0, data->X0, data->XD0); 
      vgr += sum(rmult(gradH0, data->wt0),0);
    }
    for (int i = 0; i<n; ++i) {
      gr[i] = -vgr[i]*data->parscale[i];
    }
  }

  template<class Stpm2>
  void pgrfn_SmoothLogH(int n, double * beta, double * gr, void *ex) {
    typedef SmoothLogH Smooth;
    typedef pstpm2<Smooth,Stpm2> Data;
    grfn<Data>(n, beta, gr, ex);
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec h = data->h(eta,etaD) + data->bhazard;
    mat gradh = data->gradh(eta,etaD,data->X,data->XD);
    mat Xconstraint = data->kappa * rmult(gradh,h%(h<0));
    rowvec vgr(n, fill::zeros);
    for (size_t i=0; i < data->smooth.size(); ++i) {
      SmoothLogH smoothi = data->smooth[i];
      vgr.subvec(smoothi.first_para,smoothi.last_para) += 
	(data->sp)[i] * (smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para)).t();
    }
    rowvec dconstraint = sum(Xconstraint,0);
    for (int i = 0; i<n; ++i) gr[i] += vgr[i]*data->parscale[i] + dconstraint[i];
  }

  template<class Stpm2>
  void pgrfn_SmoothHaz(int n, double * beta, double * gr, void *ex) {
    typedef SmoothHaz Smooth;
    typedef pstpm2<Smooth,Stpm2> Data;
    grfn<Data>(n, beta, gr, ex);
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    rowvec vgr(n, fill::zeros);
    for (size_t i=0; i < data->smooth.size(); ++i) {
      Smooth obj = data->smooth[i];
      vec s0 = obj.X0 * vbeta;
      vec s1 = obj.X1 * vbeta;
      vec s2 = obj.X2 * vbeta;
      vec s3 = obj.X3 * vbeta;
      vec h2 = exp(s0) % (s3 + 3*s1%s2 + s1%s1%s1);
      mat dh2sq_dbeta = 
	2*lmult(h2 % exp(s0),
		obj.X3+3*(rmult(obj.X1,s2)+rmult(obj.X2,s1))+3*lmult(s1%s1,obj.X1))+
	2*lmult(h2 % h2,obj.X0);
      vgr += (data->sp)[i]*obj.lambda*sum(lmult(obj.w, dh2sq_dbeta),0);
    }
    for (int i = 0; i<n; ++i) gr[i] += vgr[i] * data->parscale[i];
  }

  // oh, for partial template specialisation for functions...

  template<class Smooth,class Stpm2>
  double pfminfn(int n, double * beta, void *ex) { return 0.0; }

  template<>
  double pfminfn<SmoothLogH,stpm2>(int n, double * beta, void *ex) {
    return pfminfn_SmoothLogH<stpm2>(n, beta, ex); }
  template<>
  double pfminfn<SmoothLogH,stpm2PO>(int n, double * beta, void *ex) {
    return pfminfn_SmoothLogH<stpm2PO>(n, beta, ex); }
  template<>
  double pfminfn<SmoothLogH,stpm2Probit>(int n, double * beta, void *ex) {
    return pfminfn_SmoothLogH<stpm2Probit>(n, beta, ex); }
    template<>
  double pfminfn<SmoothLogH,stpm2AH>(int n, double * beta, void *ex) {
    return pfminfn_SmoothLogH<stpm2AH>(n, beta, ex); }

  template<>
  double pfminfn<SmoothHaz,stpm2>(int n, double * beta, void *ex) {
    return pfminfn_SmoothHaz<stpm2>(n, beta, ex); }
  template<>
  double pfminfn<SmoothHaz,stpm2PO>(int n, double * beta, void *ex) {
    return pfminfn_SmoothHaz<stpm2PO>(n, beta, ex); }
  template<>
  double pfminfn<SmoothHaz,stpm2Probit>(int n, double * beta, void *ex) {
    return pfminfn_SmoothHaz<stpm2Probit>(n, beta, ex); }

  template<class Smooth, class Stpm2>
  void pgrfn(int n, double * beta, double * gr, void *ex) { }

  template<>
  void pgrfn<SmoothLogH,stpm2>(int n, double * beta, double * gr, void *ex) {
    pgrfn_SmoothLogH<stpm2>(n, beta, gr, ex); }
  template<>
  void pgrfn<SmoothLogH,stpm2PO>(int n, double * beta, double * gr, void *ex) {
    pgrfn_SmoothLogH<stpm2PO>(n, beta, gr, ex); }
  template<>
  void pgrfn<SmoothLogH,stpm2Probit>(int n, double * beta, double * gr, void *ex) {
    pgrfn_SmoothLogH<stpm2Probit>(n, beta, gr, ex);}
    template<>
  void pgrfn<SmoothLogH,stpm2AH>(int n, double * beta, double * gr, void *ex) {
    pgrfn_SmoothLogH<stpm2AH>(n, beta, gr, ex); }

  template<>
  void pgrfn<SmoothHaz,stpm2>(int n, double * beta, double * gr, void *ex) {
    pgrfn_SmoothHaz<stpm2>(n, beta, gr, ex); }
  template<>
  void pgrfn<SmoothHaz,stpm2PO>(int n, double * beta, double * gr, void *ex) {
    pgrfn_SmoothHaz<stpm2PO>(n, beta, gr, ex); }
  template<>
  void pgrfn<SmoothHaz,stpm2Probit>(int n, double * beta, double * gr, void *ex) {
    pgrfn_SmoothHaz<stpm2Probit>(n, beta, gr, ex); }


  template<class Stpm2>
  SEXP optim_stpm2(SEXP args) {

    Stpm2 data(args);

    BFGS2<Stpm2> bfgs;
    bfgs.reltol = data.reltol;
    bfgs.optimWithConstraint(fminfn<Stpm2>, grfn<Stpm2>, data.init, (void *) &data, fminfn_constraint<Stpm2>);

    return List::create(_("fail")=bfgs.fail,
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));
  }

  RcppExport SEXP optim_stpm2_ph(SEXP args) {
    return optim_stpm2<stpm2>(args); }
  RcppExport SEXP optim_stpm2_po(SEXP args) {
    return optim_stpm2<stpm2PO>(args); }
  RcppExport SEXP optim_stpm2_probit(SEXP args) {
    return optim_stpm2<stpm2Probit>(args); }


  RcppExport SEXP optim_stpm2_nlm(SEXP args) {

    stpm2 data(args);

    data.init = ew_div(data.init,data.parscale);
    int n = data.init.size();
    
    Nlm nlm;
    nlm.gradtl = nlm.steptl = data.reltol;
    //nlm.method=2; nlm.stepmx=0.0;
    bool satisfied;
    do {
      nlm.optim(& fminfn_nlm<stpm2>, & grfn<stpm2>, data.init, (void *) &data);
      satisfied = fminfn_constraint<stpm2>(n,&nlm.coef[0],(void *) &data);
      if (!satisfied) data.kappa *= 2.0;
    } while (!satisfied && data.kappa<1.0e5);

    nlm.coef = ew_mult(nlm.coef, data.parscale);

    return List::create(_("itrmcd")=wrap(nlm.itrmcd),
			_("niter")=wrap(nlm.itncnt),
			_("coef")=wrap(nlm.coef),
			_("hessian")=wrap(nlm.hessian));
  }



  template<class Smooth, class Stpm2>
  double pstpm2_step_first(double logsp, void * ex) {
    typedef pstpm2<Smooth,Stpm2> Data;
    Data * data = (Data *) ex;

    data->sp[0] = exp(logsp);

    BFGS2<Data> bfgs;
    bfgs.reltol = data->reltol;

    bfgs.optimWithConstraint(pfminfn<Smooth,Stpm2>, pgrfn<Smooth,Stpm2>, data->init, ex, fminfn_constraint<Data>, false); // do not apply parscale b4/after

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, ex);

    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));
    double negll = bfgs.calc_objective(fminfn<Data>,ex);
    double gcv =  negll + data->alpha*edf;
    double bic =  negll + data->alpha*edf*log(sum(data->event));
    data->init = bfgs.coef;
    if (data->trace > 0)
      Rprintf("sp=%f\tedf=%f\tnegll=%f\tgcv=%f\tbic=%f\talpha=%f\n",data->sp[0],edf,negll,gcv,bic,data->alpha);

    return data->criterion==1 ? gcv : bic;
  }    

  template<class Smooth, class Stpm2>
  double pstpm2_step_multivariate(int n, double * logsp, void * ex) {
    typedef pstpm2<Smooth,Stpm2> Data;
    Data * data = (Data *) ex;

    double lsp = -9.0, usp = 9.0;

    for (int i=0; i < data->sp.size(); ++i)
      //data->sp[i] = exp(logsp[i]);
      data->sp[i] = exp(bound(logsp[i],lsp,usp));

    if (data->trace > 0) {
      for (int i = 0; i < data->sp.size(); ++i)
	Rprintf("sp[%i]=%f\t",i,data->sp[i]);
    }

    BFGS2<Data> bfgs; 
    bfgs.reltol = data->reltol_search;
    bfgs.optimWithConstraint(pfminfn<Smooth,Stpm2>, pgrfn<Smooth,Stpm2>, data->init, ex, fminfn_constraint<Data>, false); // do not apply parscale b4/after

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, ex);

    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));
    double negll = bfgs.calc_objective(fminfn<Data>, ex);
    double gcv =  negll + data->alpha*edf;
    double bic =  2.0*negll + data->alpha*edf*log(sum(data->event));
    data->init = bfgs.coef;

    // simple boundary constraints
    double constraint = 0.0;
    // double kappa = 1.0;
    for (int i=0; i < data->sp.size(); ++i) {
      if (logsp[i] < lsp)
	constraint +=  pow(logsp[i]-lsp,2.0);
      if (logsp[i] > usp)
	constraint +=  pow(logsp[i]-usp,2.0);
    }
    double objective =  data->criterion==1 ? 
      gcv + data->kappa/2.0*constraint : 
      bic + data->kappa/2.0*constraint;

    if (data->trace > 0) {
      Rprintf("edf=%f\tnegll=%f\tgcv=%f\tbic=%f\tobjective=%f\n",edf,negll,gcv,bic,objective);
    }
    return objective;
  }    


  template<class Smooth, class Stpm2>
  void pstpm2_step_multivariate_nlm(int n, double * logsp, double * f, void *ex) {
    *f = pstpm2_step_multivariate<Smooth,Stpm2>(n, logsp, ex);
  };


  template<class Smooth, class Stpm2>
  SEXP optim_pstpm2_first(SEXP args) { 

    typedef pstpm2<Smooth,Stpm2> Data;
    Data data(args);
    int n = data.init.size();

    for (int i=0; i<n; ++i) data.init[i] /= data.parscale[i];
    double opt_sp = exp(Brent_fmin(log(0.001),log(1000.0),&pstpm2_step_first<Smooth,Stpm2>,&data,1.0e-2));
    data.sp[0] = opt_sp;
    for (int i=0; i<n; ++i) data.init[i] *= data.parscale[i];

    BFGS2<Data> bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = data.reltol;
    bfgs.optimWithConstraint(pfminfn<Smooth,Stpm2>, pgrfn<Smooth,Stpm2>, data.init, (void *) &data, fminfn_constraint<Data>);

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, (void *) &data);
    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));

    return List::create(_("sp")=wrap(opt_sp),
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian),
          _("edf")=wrap(edf));

  }

  RcppExport SEXP optim_pstpm2LogH_first_ph(SEXP args) {
    return optim_pstpm2_first<SmoothLogH,stpm2>(args); }
  RcppExport SEXP optim_pstpm2Haz_first_ph(SEXP args) {
    return optim_pstpm2_first<SmoothHaz,stpm2>(args); }
  RcppExport SEXP optim_pstpm2LogH_first_po(SEXP args) {
    return optim_pstpm2_first<SmoothLogH,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2Haz_first_po(SEXP args) {
    return optim_pstpm2_first<SmoothHaz,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2LogH_first_probit(SEXP args) {
    return optim_pstpm2_first<SmoothLogH,stpm2Probit>(args); }
  RcppExport SEXP optim_pstpm2Haz_first_probit(SEXP args) {
    return optim_pstpm2_first<SmoothHaz,stpm2Probit>(args); }
  RcppExport SEXP optim_pstpm2LogH_first_ah(SEXP args) {
    return optim_pstpm2_first<SmoothLogH,stpm2AH>(args); }
  
  
  template<class Smooth, class Stpm2>
  SEXP optim_pstpm2_multivariate(SEXP args) {

    typedef pstpm2<Smooth,Stpm2> Data;
    Data data(args);
    int n = data.init.size();
    data.kappa = 10.0;

    NelderMead nm;
    nm.reltol = 1.0e-5;
    nm.maxit = 500;
    for (int i=0; i<n; ++i) data.init[i] /= data.parscale[i];

    NumericVector logsp(data.sp.size());
    for (int i=0; i < data.sp.size(); ++i)
      logsp[i] = log(data.sp[i]);

    bool satisfied;
    do {
      nm.optim(pstpm2_step_multivariate<Smooth,Stpm2>, logsp, (void *) &data);
      satisfied = true;
      for (int i=0; i < data.sp.size(); ++i)
	if (logsp[i]< -7.0 || logsp[i] > 7.0) satisfied = false;
      if (!satisfied) data.kappa *= 2.0;
    } while (!satisfied && data.kappa<1.0e5);

    for (int i=0; i<n; ++i) data.init[i] *= data.parscale[i];
    // for (int i=0; i < data.sp.size(); ++i)
    //   data.sp[i] = exp(logsp[i]);
    for (int i=0; i < nm.coef.size(); ++i)
      data.sp[i] = exp(nm.coef[i]);

    BFGS2<Data> bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = data.reltol;
    data.kappa = 1.0;
    bfgs.optimWithConstraint(pfminfn<Smooth,Stpm2>, pgrfn<Smooth,Stpm2>, data.init, (void *) &data, fminfn_constraint<Data>);

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, (void *) &data);
    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));

    return List::create(_("sp")=wrap(data.sp),
  			_("coef")=wrap(bfgs.coef),
  			_("hessian")=wrap(bfgs.hessian),
        _("edf")=wrap(edf));

  }

  RcppExport SEXP optim_pstpm2LogH_multivariate_ph(SEXP args) { 
    return optim_pstpm2_multivariate<SmoothLogH,stpm2>(args); }
  RcppExport SEXP optim_pstpm2Haz_multivariate_ph(SEXP args) {
    return optim_pstpm2_multivariate<SmoothHaz,stpm2>(args); }
  RcppExport SEXP optim_pstpm2LogH_multivariate_po(SEXP args) {
    return optim_pstpm2_multivariate<SmoothLogH,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2Haz_multivariate_po(SEXP args) {
    return optim_pstpm2_multivariate<SmoothHaz,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2LogH_multivariate_probit(SEXP args) {
    return optim_pstpm2_multivariate<SmoothLogH,stpm2Probit>(args); }
  RcppExport SEXP optim_pstpm2Haz_multivariate_probit(SEXP args) {
    return optim_pstpm2_multivariate<SmoothHaz,stpm2Probit>(args); }
  RcppExport SEXP optim_pstpm2LogH_multivariate_ah(SEXP args) { 
    return optim_pstpm2_multivariate<SmoothLogH,stpm2AH>(args); }

  template<class Smooth, class Stpm2>
  SEXP optim_pstpm2_multivariate_nlm(SEXP args) {

    typedef pstpm2<Smooth,Stpm2> Data;
    Data data(args);
    int n = data.init.size();
    data.kappa = 10.0;

    Nlm nlm;
    nlm.iagflg = 0;
    nlm.gradtl = 1.0e-4;
    nlm.steptl = 1.0e-4;
    nlm.msg = 9 + 4;
    for (int i=0; i<n; ++i) data.init[i] /= data.parscale[i];

    NumericVector logsp(data.sp.size());
    for (int i=0; i < data.sp.size(); ++i)
      logsp[i] = log(data.sp[i]);

    bool satisfied;
    do {
      nlm.optim(pstpm2_step_multivariate_nlm<Smooth,Stpm2>, (fcn_p) 0, logsp, (void *) &data, false);
      satisfied = true;
      for (int i=0; i < data.sp.size(); ++i)
	if (logsp[i]< -7.0 || logsp[i] > 7.0) satisfied = false;
      if (!satisfied) data.kappa *= 2.0;
    } while (!satisfied && data.kappa<1.0e5);

    for (int i=0; i<n; ++i) data.init[i] *= data.parscale[i];
    for (int i=0; i < nlm.coef.size(); ++i)
      data.sp[i] = exp(nlm.coef[i]);

    BFGS2<Data> bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = data.reltol;
    data.kappa = 1.0;
    bfgs.optimWithConstraint(pfminfn<Smooth,Stpm2>, pgrfn<Smooth,Stpm2>, data.init, (void *) &data, fminfn_constraint<Data>);

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, (void *) &data);
    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));

    return List::create(_("sp")=wrap(data.sp),
  			_("coef")=wrap(bfgs.coef),
  			_("hessian")=wrap(bfgs.hessian),
			_("edf")=wrap(edf));
  }

  RcppExport SEXP optim_pstpm2LogH_multivariate_ph_nlm(SEXP args) { 
    return optim_pstpm2_multivariate_nlm<SmoothLogH,stpm2>(args); }
  RcppExport SEXP optim_pstpm2Haz_multivariate_ph_nlm(SEXP args) {
    return optim_pstpm2_multivariate_nlm<SmoothHaz,stpm2>(args); }
  RcppExport SEXP optim_pstpm2LogH_multivariate_po_nlm(SEXP args) {
    return optim_pstpm2_multivariate_nlm<SmoothLogH,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2Haz_multivariate_po_nlm(SEXP args) {
    return optim_pstpm2_multivariate_nlm<SmoothHaz,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2LogH_multivariate_probit_nlm(SEXP args) {
    return optim_pstpm2_multivariate_nlm<SmoothLogH,stpm2Probit>(args); }
  RcppExport SEXP optim_pstpm2Haz_multivariate_probit_nlm(SEXP args) {
    return optim_pstpm2_multivariate_nlm<SmoothHaz,stpm2Probit>(args); }


  template<class Smooth, class Stpm2>
  SEXP optim_pstpm2_fixedsp(SEXP args) { 

    typedef pstpm2<Smooth,Stpm2> Data;
    Data data(args);
    BFGS2<Data> bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = data.reltol;

    bfgs.optimWithConstraint(pfminfn<Smooth,Stpm2>, pgrfn<Smooth,Stpm2>, data.init, (void *) &data, fminfn_constraint<Data>);

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, (void *) &data);
    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));

    return List::create(_("sp")=wrap(data.sp),
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian),
			_("edf")=wrap(edf));
  }

  RcppExport SEXP optim_pstpm2LogH_fixedsp_ph(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothLogH,stpm2>(args); }
  RcppExport SEXP optim_pstpm2Haz_fixedsp_ph(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothHaz,stpm2>(args); }
  RcppExport SEXP optim_pstpm2LogH_fixedsp_po(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothLogH,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2Haz_fixedsp_po(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothHaz,stpm2PO>(args); }
  RcppExport SEXP optim_pstpm2LogH_fixedsp_probit(SEXP args) { 
    return optim_pstpm2_fixedsp<SmoothLogH,stpm2Probit>(args); }
  RcppExport SEXP optim_pstpm2Haz_fixedsp_probit(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothHaz,stpm2Probit>(args); }
  RcppExport SEXP optim_pstpm2LogH_fixedsp_ah(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothLogH,stpm2AH>(args); }

  template<class Smooth, class Stpm2>
  SEXP test_pstpm2(SEXP args) { 

    pstpm2<Smooth,Stpm2> data(args);
    int n = data.init.size();

    NumericVector init = ew_div(data.init,data.parscale);

    double pfmin = pfminfn<Smooth,Stpm2>(n,&init[0],(void *) &data);
    NumericVector gr(n);
    pgrfn<Smooth,Stpm2>(n,&init[0], &gr[0], (void *) &data);

    init = ew_mult(init,data.parscale);

    return List::create(_("sp")=wrap(data.sp),
			_("pfmin")=wrap(pfmin),
			_("gr")=wrap(gr)
			);

  }

  RcppExport SEXP test_pstpm2LogH(SEXP args) {
    return test_pstpm2<SmoothLogH,stpm2>(args);
  }

  RcppExport SEXP test_pstpm2Haz(SEXP args) {
    return test_pstpm2<SmoothHaz,stpm2>(args);
  }

  // R CMD INSTALL ~/src/R/microsimulation
  // R -q -e "require(microsimulation); .Call('test_nmmin',1:2,PACKAGE='microsimulation')"

  // .Call("optim_stpm2",init,X,XD,rep(bhazard,nrow(X)),wt,ifelse(event,1,0),package="rstpm2")

} // anonymous namespace
