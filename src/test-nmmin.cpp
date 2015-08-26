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

  void Rprint(NumericVector v) {
    for (int i=0; i<v.size(); ++i) 
      Rprintf("%f ", v(i));
    Rprintf("\n");
  }

  void Rprint(vec v) {
    for (size_t i=0; i<v.size(); ++i) 
      Rprintf("%f ", v(i));
    Rprintf("\n");
  }
  void Rprint(rowvec v) {
    for (size_t i=0; i<v.size(); ++i) 
      Rprintf("%f ", v(i));
    Rprintf("\n");
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


  // object used to call class methods in the BFGS C optimiser
  // void * stpm2Object;
  // wrappers to call the stpm2Object class methods
  template<class T>
  double optimfunction(int n, double * beta, void * obj) {
    T * objT = static_cast<T*>(obj);
    NumericVector nvbeta(beta,beta+n);
    return objT->objective(nvbeta);
  }
  template<class T>
  void optimgradient(int n, double * beta, double * grad, void * obj_ptr) {
    T * obj = static_cast<T*>(obj_ptr);
    NumericVector nvbeta(beta,beta+n);
    vec gr = obj->gradient(nvbeta);
    grad = &gr[0];
  }
  template<class T>
  void optimfn_nlm(int n, double * beta, double * f, void * obj) {
    T * objT = static_cast<T*>(obj);
    NumericVector nvbeta(beta,n);
    *f = objT->objective(nvbeta);
  }

  // Various link objects
  class LogLogLink {
  public:
    vec link(vec S) { return log(-log(S)); }
    vec ilink(vec x) { return exp(-exp(x)); }
    vec h(vec eta, vec etaD) { return etaD % exp(eta); }
    vec H(vec eta, vec etaD) { return exp(eta); }
    mat gradH(vec eta, vec etaD, mat X, mat XD) { return rmult(X,exp(eta)); }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { 
      return rmult(XD, exp(eta)) + rmult(X, etaD % exp(eta));
    }
  };
  class LogitLink {
  public:
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
  class ProbitLink {
  public:
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
  class LogLink {
  public:
    vec link(vec S) { return -log(S); }
    vec ilink(vec x) { return exp(-x); }
    vec h(vec eta, vec etaD) { return etaD; }
    vec H(vec eta, vec etaD) { return eta; }
    mat gradh(vec eta, vec etaD, mat X, mat XD) { return XD; }
    mat gradH(vec eta, vec etaD, mat X, mat XD) { return X;  }
  };

  template<class Link = LogLogLink>
  class Stpm2 : public Link {
  public:
    typedef Stpm2<Link> This;
    Stpm2(SEXP sexp) {
      List list = as<List>(sexp);
      init = as<NumericVector>(list["init"]);
      X = as<mat>(list["X"]); 
      XD = as<mat>(list["XD"]); 
      bhazard = as<vec>(list["bhazard"]);
      wt = as<vec>(list["wt"]);
      event = as<vec>(list["event"]);
      time = as<vec>(list["time"]);
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
      n = init.size();
    }
    // negative log-likelihood
    double objective(NumericVector beta) {
      vec vbeta(&beta[0],n);
      vbeta = vbeta % parscale;
      vec eta = X * vbeta;
      vec etaD = XD * vbeta;
      vec h = Link::h(eta,etaD) + bhazard;
      vec H = Link::H(eta,etaD) + bhazard;
      double constraint = kappa/2.0 * sum(h % h % (h<0)); // sum(h^2 | h<0)
      vec eps = h*0.0 + 1.0e-16; 
      h = max(h,eps);
      double ll = sum(wt % event % log(h)) - sum(wt % H) - constraint;
      if (delayed == 1) {
	vec eta0 = X0 * vbeta;
	vec etaD0 = XD0 * vbeta;
	mat H0 = Link::H(eta0,etaD0);
	ll += sum(wt0 % H0);
      }
      return -ll;  
    }
    // gradient of the negative log-likelihood
    rowvec gradient(NumericVector beta) {
      vec vbeta = as<vec>(wrap(beta));
      vbeta = vbeta % parscale;
      vec eta = X * vbeta;
      vec etaD = XD * vbeta;
      vec h = Link::h(eta,etaD) + bhazard;
      // vec H = exp(eta);
      mat gradH = Link::gradH(eta,etaD,X,XD);
      mat gradh = Link::gradh(eta,etaD,X,XD);
      mat Xgrad = -gradH + rmult(gradh, event / h);
      mat Xconstraint = - kappa*rmult(gradh,h);
      vec eps = h*0.0 + 1.0e-16; // hack
      //Xgrad.rows(h<=eps) = Xconstraint.rows(h<=eps);
      Xgrad = rmult(Xgrad,wt);
      rowvec vgr = sum(Xgrad,0);
      if (delayed == 1) {
	vec eta0 = X0 * vbeta;
	vec etaD0 = XD0 * vbeta;
	mat gradH0 = Link::gradH(eta0, etaD0, X0, XD0); 
	vgr += sum(rmult(gradH0, wt0),0);
      }
      for (int i = 0; i<n; ++i) {
	// gr[i] = -vgr[i]*parscale[i];
      }
      return vgr;
    }
    bool constraint(NumericVector beta) {
      vec vbeta = as<vec>(wrap(beta));
      vbeta = vbeta % parscale;
      vec eta = X * vbeta;
      vec etaD = XD * vbeta;
      vec h = Link::h(eta, etaD) + bhazard;
      vec H = Link::H(eta, etaD);
      return all((h>0) % (H>0));
    }
    void optim(NumericVector init, bool apply_parscale = true) {
      if (apply_parscale) for (int i = 0; i<n; ++i) init[i] /= parscale[i];
      bfgs.optim(&optimfunction<This>, &optimgradient<This>, init, this);
      for (int i = 0; i<n; ++i) coef[i] *= bfgs.coef[i] * (apply_parscale ? parscale[i] : 1.0);
      bfgs.calc_hessian(&optimgradient<Stpm2<Link> >, this);
      hessian = bfgs.hessian;
    }
    void optimWithConstraint(NumericVector init, bool apply_parscale = true) {
      if (apply_parscale) for (int i = 0; i<n; ++i) init[i] /= parscale[i];
      bool satisfied;
      do {
	bfgs.optim(&optimfunction<This>, &optimgradient<This>, init, this);
	satisfied = constraint(bfgs.coef);
	if (!satisfied) kappa *= 2.0;   
      } while ((!satisfied) && kappa < 1.0e3);      
      for (int i = 0; i<n; ++i) {
	coef[i] *= bfgs.coef[i] * (apply_parscale ? parscale[i] : 1.0);
	bfgs.coef[i] = coef[i]; // to ensure that the Hessian is correctly calculated
      }
      bfgs.calc_hessian(&optimgradient<Stpm2<Link> >, this);
      hessian = as<mat>(wrap(bfgs.hessian));
    }
    NumericVector init, coef;
    mat X, XD, X0, XD0, hessian; 
    vec bhazard,wt,wt0,event,time,parscale;
    double reltol, kappa;
    int delayed, trace, n;
    BFGS bfgs;
  };

  template<class Stpm2>
  SEXP optim_stpm2(SEXP args) {
    Stpm2 model(args);
    model.optimWithConstraint(model.init);
    return List::create(_("fail")=model.bfgs.fail, 
			_("coef")=wrap(model.coef),
			_("hessian")=wrap(model.hessian));
  }
  RcppExport SEXP optim_stpm2_ph(SEXP args) {
    return optim_stpm2<Stpm2<LogLogLink> >(args); }
  RcppExport SEXP optim_stpm2_po(SEXP args) {
    return optim_stpm2<Stpm2<LogitLink> >(args); }
  RcppExport SEXP optim_stpm2_probit(SEXP args) {
    return optim_stpm2<Stpm2<ProbitLink> >(args); }
  RcppExport SEXP optim_stpm2_ah(SEXP args) {
    return optim_stpm2<Stpm2<LogLink> >(args); }
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

  // penalised smoothers (currently _not_ sub-classes of Pstpm2)
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
    double penalty(NumericVector beta, vec parscale, vec sp) {
      double value = 0.0;
      int n = beta.size();
      vec vbeta(&beta[0],n);
      vbeta = vbeta % parscale;
      for (size_t i=0; i < smooth.size(); ++i) {
	Smoother smoothi = smooth[i];
	value += sp[i]/2 * 
	  dot(vbeta.subvec(smoothi.first_para,smoothi.last_para),
	      smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para));
      }
      return value;  
    }
    vec penalty_gradient(NumericVector beta, vec parscale, NumericVector sp) {
      int n = beta.size();
      vec vbeta(&beta[0],n);
      vbeta = vbeta % parscale;
      rowvec vgr(n, fill::zeros);
      for (size_t i=0; i < smooth.size(); ++i) {
	Smoother smoothi = smooth[i];
	vgr.subvec(smoothi.first_para,smoothi.last_para) += 
	  sp[i] * (smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para)).t();
      }
      for (int i = 0; i<n; ++i) vgr[i] += vgr[i]*parscale[i];
      return vgr;
    }
    std::vector<Smoother> smooth;
  };
  class SmoothHaz {
  public:
    struct Smoother {
      mat X0, X1, X2, X3;
      vec w;
      double lambda;
    };
    SmoothHaz(SEXP sexp) {
      List list = as<List>(sexp);
      List lsmooth = as<List>(list["smooth"]);
      for(int i=0; i<lsmooth.size(); ++i) {
	List lsmoothi = as<List>(lsmooth[i]);
	Smoother smoothi =  {
	  as<mat>(lsmoothi["X0"]), 
	  as<mat>(lsmoothi["X1"]), 
	  as<mat>(lsmoothi["X2"]), 
	  as<mat>(lsmoothi["X3"]), 
	  as<vec>(lsmoothi["w"]), 
	  as<double>(lsmoothi["lambda"])
	};
	smooth.push_back(smoothi);
      }
    }
    double penalty(NumericVector beta, vec parscale, NumericVector sp) {
      double value = 0.0;
      int n = beta.size();
      vec vbeta(&beta[0],n);
      vbeta = vbeta % parscale;
      for (size_t i=0; i < smooth.size(); ++i) {
	Smoother obj = smooth[i];
	vec s0 = obj.X0 * vbeta;
	vec s1 = obj.X1 * vbeta;
	vec s2 = obj.X2 * vbeta;
	vec s3 = obj.X3 * vbeta;
	vec h2 = exp(s0) % (s3 + 3*s1%s2 + s1%s1%s1);
	value += sp[i]/2 * obj.lambda * sum(obj.w % h2 % h2);
      }
      return value;
    }
  vec penalty_gradient(NumericVector beta, vec parscale, NumericVector sp) {
    int n = beta.size();
    vec vbeta(&beta[0],n);
    vbeta = vbeta % parscale;
    rowvec vgr(n, fill::zeros);
    for (size_t i=0; i < smooth.size(); ++i) {
      Smoother obj = smooth[i];
      vec s0 = obj.X0 * vbeta;
      vec s1 = obj.X1 * vbeta;
      vec s2 = obj.X2 * vbeta;
      vec s3 = obj.X3 * vbeta;
      vec h2 = exp(s0) % (s3 + 3*s1%s2 + s1%s1%s1);
      mat dh2sq_dbeta = 
	2*lmult(h2 % exp(s0),
		obj.X3+3*(rmult(obj.X1,s2)+rmult(obj.X2,s1))+3*lmult(s1%s1,obj.X1))+
	2*lmult(h2 % h2,obj.X0);
      vgr += sp[i]*obj.lambda*sum(lmult(obj.w, dh2sq_dbeta),0);
    }
    for (int i = 0; i<n; ++i) vgr[i] += vgr[i] * parscale[i];
    return vgr;
  }
    std::vector<Smoother> smooth;
  };

  template<class T>
  double pstpm2_multivariate_step(int n, double * logsp_ptr, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    NumericVector logsp(logsp_ptr,logsp_ptr+n);
    return model->multivariate_step(logsp);
  }    
  template<class T>
  double pstpm2_first_step(double logsp, void * model_ptr) {
    T * model = static_cast<T *>(model_ptr);
    return model->first_step(logsp);
  }    

  // Penalised link-based models
  template<class Stpm2Type = Stpm2<LogLogLink>, class Smooth = SmoothLogH>
  class Pstpm2 : public Stpm2Type, public Smooth {
  public:
    typedef Pstpm2<Stpm2Type> Pstpm2Type;
    Pstpm2(SEXP sexp) : Stpm2Type(sexp), Smooth(sexp) {
      List list = as<List>(sexp);
      sp = as<NumericVector>(list["sp"]);
      reltol_search = as<double>(list["reltol_search"]);
      alpha = as<double>(list["alpha"]);
      criterion = as<int>(list["criterion"]);
      parscale = as<vec>(list["parscale"]);
    }
    double objective(NumericVector beta) {
      return this->objective(beta) - Smooth::penalty(beta,parscale,sp);
    }
    vec gradient(NumericVector beta) {
      return this->gradient(beta) - Smooth::penalty_gradient(beta,parscale,sp);
    }
    double first_step(double logsp) {
      sp[0] = exp(logsp);
      this->optimWithConstraint(this->init, false); // do not apply parscale b4/after
      NumericMatrix hessian0 = this->bfgs.calc_hessian(&optimgradient<Stpm2Type>, (void *) this);
      if (this->trace > 0)  {
	Rprintf("Debug on trace calculation. Coef:\n");
	Rprint(this->bfgs.coef);
	Rprintf("Hessian0:\n");
	Rprint(hessian0);
	Rprintf("Hessian:\n");
	Rprint(this->bfgs.hessian);
      }
      double edf = arma::trace(solve(as<mat>(this->bfgs.hessian),as<mat>(hessian0)));
      double negll = this->bfgs.calc_objective(&optimfunction<Stpm2Type>, (void *) this);
      double gcv =  negll + alpha*edf;
      double bic =  negll + alpha*edf*log(sum(this->event));
      this->init = this->bfgs.coef;
      if (this->trace > 0)
	Rprintf("sp=%f\tedf=%f\tnegll=%f\tgcv=%f\tbic=%f\talpha=%f\n",sp[0],edf,negll,gcv,bic,alpha);
      return criterion==1 ? gcv : bic;
    }
    double multivariate_step(NumericVector logsp) {
      double lsp = -9.0, usp = 9.0;
      for (int i=0; i < sp.size(); ++i)
	sp[i] = exp(bound(logsp[i],lsp,usp));
      if (this->trace > 0) {
	for (int i = 0; i < sp.size(); ++i)
	  Rprintf("sp[%i]=%f\t",i,sp[i]);
      }
      this->optimWithConstraint(this->init, false); // do not apply parscale b4/after
      NumericMatrix hessian0 = this->bfgs.calc_hessian(&optimgradient<Stpm2Type>, (void *) this);
      double edf = arma::trace(solve(as<mat>(this->bfgs.hessian),as<mat>(hessian0)));
      double negll = this->bfgs.calc_objective(&optimfunction<Stpm2Type>, (void *) this);
      double gcv =  negll + alpha*edf;
      double bic =  2.0*negll + alpha*edf*log(sum(this->event));
      this->init = this->bfgs.coef;
      // simple boundary constraints
      double constraint = 0.0;
      // double kappa = 1.0;
      for (int i=0; i < sp.size(); ++i) {
	if (logsp[i] < lsp)
	  constraint +=  pow(logsp[i]-lsp,2.0);
	if (logsp[i] > usp)
	  constraint +=  pow(logsp[i]-usp,2.0);
      }
      double objective =  criterion==1 ? 
	gcv + this->kappa/2.0*constraint : 
	bic + this->kappa/2.0*constraint;
      if (this->trace > 0) {
	Rprintf("edf=%f\tnegll=%f\tgcv=%f\tbic=%f\tobjective=%f\n",edf,negll,gcv,bic,objective);
      }
      return objective;
    }
    SEXP optim_fixed() { 
      // for (int i=0; i<n; ++i) init[i] *= parscale[i];
      this->optimWithConstraint(this->init);
      NumericMatrix hessian0 = this->bfgs.calc_hessian(&optimgradient<Stpm2Type>, (void *) this);
      double edf = trace(solve(as<mat>(this->bfgs.hessian),as<mat>(hessian0)));
      return List::create(_("sp")=wrap(sp),
			  _("coef")=wrap(this->bfgs.coef),
			  _("hessian")=wrap(this->bfgs.hessian),
			  _("edf")=wrap(edf));
    }
    SEXP optim_first() { 
      for (int i=0; i < this->n; ++i) this->init[i] /= this->parscale[i];
      double opt_sp = exp(Brent_fmin(log(0.001),log(1000.0),&(pstpm2_first_step<Pstpm2Type>),(void *) this,1.0e-2));
      sp[0] = opt_sp;
      for (int i=0; i < this->n; ++i) this->init[i] *= this->parscale[i];
      return optim_fixed();
    }
    SEXP optim_multivariate() {
      this->kappa = 10.0;
      NelderMead nm;
      nm.reltol = 1.0e-5;
      nm.maxit = 500;
      nm.hessianp = false;
      for (int i=0; i < this->n; ++i) this->init[i] /= this->parscale[i];
      NumericVector logsp(sp.size());
      for (int i=0; i < sp.size(); ++i)
	logsp[i] = log(sp[i]);
      bool satisfied;
      do {
	nm.optim(&pstpm2_multivariate_step<Pstpm2Type>, logsp, (void *) this);
	satisfied = true;
	for (int i=0; i < sp.size(); ++i)
	  if (logsp[i]< -7.0 || logsp[i] > 7.0) satisfied = false;
	if (!satisfied) this->kappa *= 2.0;
      } while (!satisfied && this->kappa<1.0e5);
      for (int i=0; i < this->n; ++i) this->init[i] *= this->parscale[i];
      for (int i=0; i < nm.coef.size(); ++i)
	sp[i] = exp(nm.coef[i]);
      this->bfgs.coef = this->init;     // ?
      // this->bfgs.reltol = reltol; // ?
      return optim_fixed();
    }
    NumericVector sp;
    vec parscale;
    double alpha, reltol_search; 
    int criterion;
  };

  // template<class Smooth, class Stpm2>
  // void pstpm2_step_multivariate_nlm(int n, double * logsp, double * f, void *ex) {
  //   *f = pstpm2_step_multivariate<Smooth,Stpm2>(n, logsp, ex);
  // };

  RcppExport SEXP optim_pstpm2LogH_first_ph(SEXP args) {
    Pstpm2<Stpm2<LogLogLink>,SmoothLogH> model(args);
    return model.optim_first();
  }
  RcppExport SEXP optim_pstpm2Haz_first_ph(SEXP args) {
    Pstpm2<Stpm2<LogLogLink>,SmoothHaz> model(args);
    return model.optim_first();
  }
  RcppExport SEXP optim_pstpm2LogH_first_po(SEXP args) {
    Pstpm2<Stpm2<LogitLink>,SmoothLogH> model(args);
    return model.optim_first();
  }
  RcppExport SEXP optim_pstpm2Haz_first_po(SEXP args) {
    Pstpm2<Stpm2<LogitLink>,SmoothHaz> model(args);
    return model.optim_first();
  }
  RcppExport SEXP optim_pstpm2LogH_first_probit(SEXP args) {
    Pstpm2<Stpm2<ProbitLink>,SmoothLogH> model(args);
    return model.optim_first();
  }
  RcppExport SEXP optim_pstpm2Haz_first_probit(SEXP args) {
    Pstpm2<Stpm2<ProbitLink>,SmoothHaz> model(args);
    return model.optim_first();
  }
  RcppExport SEXP optim_pstpm2LogH_first_ah(SEXP args) {
    Pstpm2<Stpm2<LogLink>,SmoothLogH> model(args);
    return model.optim_first();
  }

  RcppExport SEXP optim_pstpm2LogH_multivariate_ph(SEXP args) {
    Pstpm2<Stpm2<LogLogLink>,SmoothLogH> model(args);
    return model.optim_multivariate();
  }
  RcppExport SEXP optim_pstpm2Haz_multivariate_ph(SEXP args) {
    Pstpm2<Stpm2<LogLogLink>,SmoothHaz> model(args);
    return model.optim_multivariate();
  }
  RcppExport SEXP optim_pstpm2LogH_multivariate_po(SEXP args) {
    Pstpm2<Stpm2<LogitLink>,SmoothLogH> model(args);
    return model.optim_multivariate();
  }
  RcppExport SEXP optim_pstpm2Haz_multivariate_po(SEXP args) {
    Pstpm2<Stpm2<LogitLink>,SmoothHaz> model(args);
    return model.optim_multivariate();
  }
  RcppExport SEXP optim_pstpm2LogH_multivariate_probit(SEXP args) {
    Pstpm2<Stpm2<ProbitLink>,SmoothLogH> model(args);
    return model.optim_multivariate();
  }
  RcppExport SEXP optim_pstpm2Haz_multivariate_probit(SEXP args) {
    Pstpm2<Stpm2<ProbitLink>,SmoothHaz> model(args);
    return model.optim_multivariate();
  }
  RcppExport SEXP optim_pstpm2LogH_multivariate_ah(SEXP args) {
    Pstpm2<Stpm2<LogLink>,SmoothLogH> model(args);
    return model.optim_multivariate();
  }



  // template<class Smooth, class Stpm2>
  // SEXP optim_pstpm2_fixedsp(SEXP args) { 

  //   typedef pstpm2<Smooth,Stpm2> Data;
  //   Data data(args);
  //   BFGS2<Data> bfgs;
  //   bfgs.coef = data.init;
  //   bfgs.reltol = data.reltol;

  //   bfgs.optimWithConstraint(pfminfn<Smooth,Stpm2>, pgrfn<Smooth,Stpm2>, data.init, (void *) &data, fminfn_constraint<Data>);

  //   NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, (void *) &data);
  //   double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));

  //   return List::create(_("sp")=wrap(data.sp),
  // 			_("coef")=wrap(bfgs.coef),
  // 			_("hessian")=wrap(bfgs.hessian),
  // 			_("edf")=wrap(edf));
  // }

  // RcppExport SEXP optim_pstpm2LogH_fixedsp_ph(SEXP args) {
  //   return optim_pstpm2_fixedsp<SmoothLogH,stpm2>(args); }
  // RcppExport SEXP optim_pstpm2Haz_fixedsp_ph(SEXP args) {
  //   return optim_pstpm2_fixedsp<SmoothHaz,stpm2>(args); }
  // RcppExport SEXP optim_pstpm2LogH_fixedsp_po(SEXP args) {
  //   return optim_pstpm2_fixedsp<SmoothLogH,stpm2PO>(args); }
  // RcppExport SEXP optim_pstpm2Haz_fixedsp_po(SEXP args) {
  //   return optim_pstpm2_fixedsp<SmoothHaz,stpm2PO>(args); }
  // RcppExport SEXP optim_pstpm2LogH_fixedsp_probit(SEXP args) { 
  //   return optim_pstpm2_fixedsp<SmoothLogH,stpm2Probit>(args); }
  // RcppExport SEXP optim_pstpm2Haz_fixedsp_probit(SEXP args) {
  //   return optim_pstpm2_fixedsp<SmoothHaz,stpm2Probit>(args); }
  // RcppExport SEXP optim_pstpm2LogH_fixedsp_ah(SEXP args) {
  //   return optim_pstpm2_fixedsp<SmoothLogH,stpm2AH>(args); }

  // template<class Smooth, class Stpm2>
  // SEXP test_pstpm2(SEXP args) { 

  //   pstpm2<Smooth,Stpm2> data(args);
  //   int n = data.init.size();

  //   NumericVector init = ew_div(data.init,data.parscale);

  //   double pfmin = pfminfn<Smooth,Stpm2>(n,&init[0],(void *) &data);
  //   NumericVector gr(n);
  //   pgrfn<Smooth,Stpm2>(n,&init[0], &gr[0], (void *) &data);

  //   init = ew_mult(init,data.parscale);

  //   return List::create(_("sp")=wrap(data.sp),
  // 			_("pfmin")=wrap(pfmin),
  // 			_("gr")=wrap(gr)
  // 			);

  // }

  // RcppExport SEXP test_pstpm2LogH(SEXP args) {
  //   return test_pstpm2<SmoothLogH,stpm2>(args);
  // }

  // RcppExport SEXP test_pstpm2Haz(SEXP args) {
  //   return test_pstpm2<SmoothHaz,stpm2>(args);
  // }


  /** 
      Extension of stpm2 and pstpm2 to include shared frailties with a cluster variable
  **/
  template<class Stpm2Type>
  class Stpm2Shared : public Stpm2Type {
  public:
    typedef std::vector<int> Index;
    typedef std::map<int,Index> IndexMap;
    Stpm2Shared(SEXP sexp) : Stpm2Type(sexp) {
      List list = as<List>(sexp);
      IntegerVector cluster = as<IntegerVector>(list["cluster"]);
      // wragged array indexed by a map of vectors
      for (int i=0; i<cluster.size(); ++i) {
	clusters[cluster[i]].push_back(i);
      }
    }
    /** 
	Objective function for a Gamma shared frailty
	Assumes that weights == 1.
	The implementation is incomplete, requiring either the use of optimisation without
	gradients or the calculation of the gradients.
    **/
    double objective(NumericVector beta) {
      int n = beta.size();
      vec vbeta(&beta[0],n-1); // theta is the last parameter in beta
      double theta = beta[n-1];
      vbeta = vbeta % this->parscale;
      vec eta = this->X * vbeta;
      vec etaD = this->XD * vbeta;
      vec h = this->h(eta,etaD) + this->bhazard;
      vec H = this->H(eta,etaD) + this->bhazard;
      double constraint = this->kappa/2.0 * sum(h % h % (h<0)); // sum(h^2 | h<0)
      vec eps = h*0.0 + 1.0e-16; 
      h = max(h,eps);
      mat H0;
      //double ll = sum(wt % event % log(h)) - sum(wt % H) - constraint;
      if (this->delayed == 1) {
	vec eta0 = this->X0 * vbeta;
	vec etaD0 = this->XD0 * vbeta;
	H0 = this->H(eta0,etaD0);
	//ll += sum(wt0 % H0);
      }
      double ll;
      for (IndexMap::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
	int mi=0;
	double sumH = 0.0, sumHenter = 0.0;
	for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	  if (this->event(*j)==1) {
	    mi++;
	    ll += log(this->h(*j));
	  }
	  sumH += this->H(*j);
	  if (this->delayed == 1)
	    sumHenter += this->H0(*j);
	}
	ll -= (1.0/theta+mi)*log(1.0+theta*sumH);
	if (this->delayed == 1)
	  ll += 1.0/theta*log(1.0+theta*sumHenter);
	if (mi>0) {
	  int k=1;
	  for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	    if (this->event(*j)==1) {
	      ll += log(1.0+theta*(mi-k));
	      k++;
	    }
	  }
	}
      }
      ll -= constraint;
      return -ll;  
    }
    IndexMap clusters;
  };


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

  /** 
      Utility function for calculating a gamma frailty log-likelihood with delayed entry within R code 
   **/
  RcppExport SEXP llgammafrailtydelayed(SEXP _theta, SEXP _h, SEXP _H, SEXP _H0, SEXP _d,
				    SEXP _cluster) {
    double theta = as<double>(_theta);
    NumericVector h = as<NumericVector>(_h);
    NumericVector H = as<NumericVector>(_H);
    NumericVector H0 = as<NumericVector>(_H0); // assumes that _H0 is zero for data not left truncated
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
      double sumH = 0.0, sumHenter = 0.0;
      for (Index::iterator j=it->second.begin(); j!=it->second.end(); ++j) {
	if (d[*j]==1) {
	  mi++;
	  ll += log(h[*j]);
	}
	sumH += H[*j];
	sumHenter += H0[*j];
      }
      ll -= (1.0/theta+mi)*log(1.0+theta*sumH);
      ll += 1.0/theta*log(1.0+theta*sumHenter);
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

  // .Call("optim_stpm2",init,X,XD,rep(bhazard,nrow(X)),wt,ifelse(event,1,0),package="rstpm2")

} // anonymous namespace
