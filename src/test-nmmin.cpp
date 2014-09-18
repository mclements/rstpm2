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
    return as<vec>(wrap(pnorm(as<NumericVector>(wrap<vec>(x)),0.0,1.0)));
  }

  vec dnorm01(vec x) {
    return as<vec>(wrap(dnorm(as<NumericVector>(wrap<vec>(x)),0.0,1.0)));
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
    mat X, XD,X0,XD0; 
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

  // class stpm2Probit : public stpm2 {
  // public:
  //   stpm2Probit(SEXP sexp) : stpm2(sexp) { }
  //   vec link(vec S) { return -logit(S); }
  //   vec ilink(vec x) { return expit(-x); }
  //   vec h(vec eta, vec etaD) { return etaD % exp(eta) % expit(-eta); }
  //   vec H(vec eta, vec etaD) { return log(1+exp(eta)); }
  //   mat gradH(vec eta, vec etaD, mat X, mat XD) { 
  //     return rmult(X,exp(eta) % expit(-eta));
  //   }
  //   mat gradh(vec eta, vec etaD, mat X, mat XD) { 
  //     return rmult(X, etaD % exp(eta) % expit(-eta)) - 
  // 	rmult(X, exp(2*eta) % etaD % expit(-eta) % expit(-eta)) +
  // 	rmult(XD, exp(eta) % expit(-eta));
  //   }
  // };

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

  template<class Smooth>
  class pstpm2 : public stpm2 {
  public:
    pstpm2(SEXP sexp) : stpm2(sexp) {
      List list = as<List>(sexp);
      List lsmooth = as<List>(list["smooth"]);
      sp = as<NumericVector>(list["sp"]);
      reltol_search = as<double>(list["reltol_search"]);
      reltol = as<double>(list["reltol"]);
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
      } while ((!satisfied) && data->kappa < 1.0e5);
      if (apply_parscale) for (int i = 0; i<n; ++i) coef[i] *= data->parscale[i];
      hessian = calc_hessian(gr, ex, eps);
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
      eta = data->X0 * vbeta;
      etaD = data->XD0 * vbeta;
      H = data->H(eta,etaD);
      ll += sum(data->wt0 % H);
    }
    return -ll;  
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


  // empty default
  template<class Smooth>
  double pfminfn(int n, double * beta, void *ex) { return 0.0; };

  // log H penalty - using penalty matrices from mgcv
  template<>
  double pfminfn<SmoothLogH>(int n, double * beta, void *ex) {
    typedef SmoothLogH Smooth;
    typedef pstpm2<Smooth> Data;
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
  template<>
  double pfminfn<SmoothHaz>(int n, double * beta, void *ex) {
    typedef SmoothHaz Smooth;
    typedef pstpm2<Smooth> Data;
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
      mat gradH0 = data->gradH(data->X0 * vbeta, data->XD0 * vbeta, data->X0, data->XD0); 
      vgr += sum(rmult(gradH0, data->wt0),0);
    }
    for (int i = 0; i<n; ++i) {
      gr[i] = -vgr[i]*data->parscale[i];
    }
  }

  template<class Smooth>
  void pgrfn(int n, double * beta, double * gr, void *ex) { };

  template<>
  void pgrfn<SmoothLogH>(int n, double * beta, double * gr, void *ex) {
    typedef SmoothLogH Smooth;
    typedef pstpm2<Smooth> Data;
    grfn<Data>(n, beta, gr, ex);
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    rowvec vgr(n, fill::zeros);
    for (size_t i=0; i < data->smooth.size(); ++i) {
      SmoothLogH smoothi = data->smooth[i];
      vgr.subvec(smoothi.first_para,smoothi.last_para) += 
	(data->sp)[i] * (smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para)).t();
    }
    for (int i = 0; i<n; ++i) gr[i] += vgr[i]*data->parscale[i];
  }

  template<>
  void pgrfn<SmoothHaz>(int n, double * beta, double * gr, void *ex) {
    typedef SmoothHaz Smooth;
    typedef pstpm2<Smooth> Data;
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


  RcppExport SEXP optim_stpm2(SEXP args) {

    stpm2 data(args);

    BFGS2<stpm2> bfgs;
    bfgs.reltol = data.reltol;
    bfgs.optimWithConstraint(fminfn<stpm2>, grfn<stpm2>, data.init, (void *) &data, fminfn_constraint<stpm2>);

    return List::create(_("fail")=bfgs.fail,
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));
  }


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


  template<class Smooth>
  double pstpm2_step_first(double logsp, void * ex) {
    typedef pstpm2<Smooth> Data;
    Data * data = (Data *) ex;

    data->sp[0] = exp(logsp);

    BFGS2<Data> bfgs;
    bfgs.reltol = data->reltol;

    bfgs.optimWithConstraint(pfminfn<Smooth>, pgrfn<Smooth>, data->init, ex, fminfn_constraint<Data>, false); // do not apply parscale b4/after

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

  template<class Smooth>
  double pstpm2_step_multivariate(int n, double * logsp, void * ex) {
    typedef pstpm2<Smooth> Data;
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
    bfgs.optimWithConstraint(pfminfn<Smooth>, pgrfn<Smooth>, data->init, ex, fminfn_constraint<Data>, false); // do not apply parscale b4/after

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


  template<class Smooth>
  void pstpm2_step_multivariate_nlm(int n, double * logsp, double * f, void *ex) {
    *f = pstpm2_step_multivariate<Smooth>(n, logsp, ex);
  };


  template<class Smooth>
  SEXP optim_pstpm2_first(SEXP args) { 

    typedef pstpm2<Smooth> Data;
    Data data(args);
    int n = data.init.size();

    for (int i=0; i<n; ++i) data.init[i] /= data.parscale[i];
    double opt_sp = exp(Brent_fmin(log(0.001),log(1000.0),&pstpm2_step_first<Smooth>,&data,1.0e-2));
    data.sp[0] = opt_sp;
    for (int i=0; i<n; ++i) data.init[i] *= data.parscale[i];

    BFGS2<Data> bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = data.reltol;
    bfgs.optimWithConstraint(pfminfn<Smooth>, pgrfn<Smooth>, data.init, (void *) &data, fminfn_constraint<Data>);

    return List::create(_("sp")=wrap(opt_sp),
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));

  }

  RcppExport SEXP optim_pstpm2LogH_first(SEXP args) {
    return optim_pstpm2_first<SmoothLogH>(args);
  }

  RcppExport SEXP optim_pstpm2Haz_first(SEXP args) {
    return optim_pstpm2_first<SmoothHaz>(args);
  }

  template<class Smooth>
  SEXP optim_pstpm2_multivariate(SEXP args) {

    typedef pstpm2<Smooth> Data;
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
      nm.optim(pstpm2_step_multivariate<Smooth>, logsp, (void *) &data);
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
    bfgs.optimWithConstraint(pfminfn<Smooth>, pgrfn<Smooth>, data.init, (void *) &data, fminfn_constraint<Data>);

    return List::create(_("sp")=wrap(data.sp),
  			_("coef")=wrap(bfgs.coef),
  			_("hessian")=wrap(bfgs.hessian));

  }

  template<class Smooth>
  SEXP optim_pstpm2_multivariate_nlm(SEXP args) {

    typedef pstpm2<Smooth> Data;
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
      nlm.optim(pstpm2_step_multivariate_nlm<Smooth>, (fcn_p) 0, logsp, (void *) &data, false);
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
    bfgs.optimWithConstraint(pfminfn<Smooth>, pgrfn<Smooth>, data.init, (void *) &data, fminfn_constraint<Data>);

    return List::create(_("sp")=wrap(data.sp),
  			_("coef")=wrap(bfgs.coef),
  			_("hessian")=wrap(bfgs.hessian));
  }


  RcppExport SEXP optim_pstpm2LogH_multivariate(SEXP args) {
    return optim_pstpm2_multivariate_nlm<SmoothLogH>(args);
  }

  RcppExport SEXP optim_pstpm2Haz_multivariate(SEXP args) {
    return optim_pstpm2_multivariate_nlm<SmoothHaz>(args);
  }



  template<class Smooth>
  SEXP optim_pstpm2_fixedsp(SEXP args) { 

    typedef pstpm2<Smooth> Data;
    Data data(args);
    BFGS2<Data> bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = data.reltol;

    bfgs.optimWithConstraint(pfminfn<Smooth>, pgrfn<Smooth>, data.init, (void *) &data, fminfn_constraint<Data>);

    return List::create(_("sp")=wrap(data.sp),
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));

  }

  RcppExport SEXP optim_pstpm2LogH_fixedsp(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothLogH>(args);
  }

  RcppExport SEXP optim_pstpm2Haz_fixedsp(SEXP args) {
    return optim_pstpm2_fixedsp<SmoothHaz>(args);
  }



  template<class Smooth>
  SEXP test_pstpm2(SEXP args) { 

    pstpm2<Smooth> data(args);
    int n = data.init.size();

    NumericVector init = ew_div(data.init,data.parscale);

    double pfmin = pfminfn<Smooth>(n,&init[0],(void *) &data);
    NumericVector gr(n);
    pgrfn<Smooth>(n,&init[0], &gr[0], (void *) &data);

    init = ew_mult(init,data.parscale);

    return List::create(_("sp")=wrap(data.sp),
			_("pfmin")=wrap(pfmin),
			_("gr")=wrap(gr)
			);

  }

  RcppExport SEXP test_pstpm2LogH(SEXP args) {
    return test_pstpm2<SmoothLogH>(args);
  }

  RcppExport SEXP test_pstpm2Haz(SEXP args) {
    return test_pstpm2<SmoothHaz>(args);
  }


  /* PROPORTIONAL ODDS MODELS */

  RcppExport SEXP optim_stpm2_po(SEXP args) {

    stpm2PO data(args);

    BFGS2<stpm2PO> bfgs;
    bfgs.reltol = data.reltol;
    bfgs.optimWithConstraint(fminfn<stpm2PO>, grfn<stpm2PO>, data.init, (void *) &data, fminfn_constraint<stpm2PO>);

    return List::create(_("fail")=bfgs.fail,
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));
  }


  /* PROBIT MODELS */

  template<class Data>
  double fminfn_probit(int n, double * beta, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec H = -log(pnorm01(-eta));
    vec h =  dnorm01(eta) / pnorm01(-eta) % etaD + data->bhazard;
    double constraint = data->kappa/2.0 * sum(h % h % (h<0)); // sum(h^2 | h<0)
    vec eps = h*0.0 + 1.0e-16; 
    h = max(h,eps);
    double ll = sum(data->wt % data->event % log(h)) - sum(data->wt % H) - constraint;
    if (data->delayed == 1) {
      vec eta0 = data->X0 * vbeta;
      vec H0 = -log(pnorm01(-eta0));
      ll += sum(data->wt0 % H0);
    }
    return -ll;  
  }

  template<class Data>
  bool fminfn_probit_constraint(int n, double * beta, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec h =  dnorm01(eta) / pnorm01(-eta) % etaD + data->bhazard;
    return all(h>0);
  }

  template<class Data>
  void fminfn_probit_nlm(int n, double * beta, double * f, void *ex) {
    *f = fminfn_po<Data>(n, beta, ex);
  };

  template<class Data>
  void grfn_probit(int n, double * beta, double * gr, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vbeta = vbeta % data->parscale;
    vec eta = data->X * vbeta;
    vec etaD = data->XD * vbeta;
    vec H = -log(pnorm01(-eta));
    vec h =  dnorm01(eta) / pnorm01(-eta) % etaD + data->bhazard;
    mat gradH = rmult(data->X, dnorm01(eta)/pnorm01(-eta)); 
    mat gradh = rmult(data->X, -eta % dnorm01(eta) % etaD / pnorm01(-eta)) +
      rmult(data->X, dnorm01(eta) % dnorm01(eta)/pnorm01(-eta)/pnorm01(-eta) % etaD) +
      rmult(data->XD, dnorm01(eta) / pnorm01(-eta));
    mat Xgrad = -gradH + rmult(gradh, data->event/h);
    mat Xconstraint = rmult(gradh, data->kappa * h);
    vec eps = h*0.0 + 1.0e-16; // hack
    Xgrad.rows(h<=eps) = Xconstraint.rows(h<=eps);
    Xgrad = rmult(Xgrad,data->wt);
    rowvec vgr = sum(Xgrad,0);
    if (data->delayed == 1) {
      vec eta0 = data->X0 * vbeta;
      mat gradH0 = rmult(data->X, data->wt0 % dnorm01(eta0)/pnorm01(-eta0)); 
      vgr += sum(-gradH0,0);
    }
    for (int i = 0; i<n; ++i) {
      gr[i] = -vgr[i]*data->parscale[i];
    }
  }

  RcppExport SEXP optim_stpm2_probit(SEXP args) {

    stpm2 data(args);

    BFGS2<stpm2> bfgs;
    bfgs.reltol = data.reltol;
    bfgs.optimWithConstraint(fminfn_probit<stpm2>, grfn_probit<stpm2>, data.init, (void *) &data, fminfn_probit_constraint<stpm2>);

    return List::create(_("fail")=bfgs.fail,
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));
  }


  RcppExport SEXP optim_stpm2_probit_nlm(SEXP args) {

    stpm2 data(args);

    data.init = ew_div(data.init,data.parscale);
    int n = data.init.size();
    
    Nlm nlm;
    nlm.gradtl = nlm.steptl = data.reltol;
    //nlm.method=2; nlm.stepmx=0.0;
    bool satisfied;
    do {
      nlm.optim(& fminfn_probit_nlm<stpm2>, & grfn_probit<stpm2>, data.init, (void *) &data);
      satisfied = fminfn_probit_constraint<stpm2>(n,&nlm.coef[0],(void *) &data);
      if (!satisfied) data.kappa *= 2.0;
    } while (!satisfied && data.kappa<1.0e5);

    nlm.coef = ew_mult(nlm.coef, data.parscale);

    return List::create(_("itrmcd")=wrap(nlm.itrmcd),
			_("niter")=wrap(nlm.itncnt),
			_("coef")=wrap(nlm.coef),
			_("hessian")=wrap(nlm.hessian));
  }





  // R CMD INSTALL ~/src/R/microsimulation
  // R -q -e "require(microsimulation); .Call('test_nmmin',1:2,PACKAGE='microsimulation')"

  // .Call("optim_stpm2",init,X,XD,rep(bhazard,nrow(X)),wt,ifelse(event,1,0),package="rstpm2")

} // anonymous namespace
