#include <R_ext/Applic.h>
#include <RcppArmadillo.h>
#include <vector>
#include <map>
#include <float.h> /* DBL_EPSILON */

namespace {

  using namespace Rcpp;
  using namespace arma;

  typedef double optimfn(int, double *, void *);
  typedef void optimgr(int, double *, double *, void *);

  double min(double a, double b) { return a < b ? a : b; }
  double max(double a, double b) { return a < b ? b : a; }
  double bound(double x, double lower, double upper) { return x < lower ? lower : (x > upper ? upper : x); }

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

  // void nmmin(int n, double *Bvec, double *X, double *Fmin, optimfn fminfn,
  // 	   int *fail, double abstol, double intol, void *ex,
  // 	   double alpha, double bet, double gamm, int trace,
  // 	   int *fncount, int maxit)

  class NelderMead {
  public:
    NelderMead(int trace = 0, int maxit = 500, 
	       double abstol = - INFINITY,
	       double reltol = 1.0e-8, 
	       double alpha = 1.0, double beta = 0.5, double gamma = 2.0) : 
      trace(trace), maxit(maxit), abstol(abstol), reltol(reltol), 
      alpha(alpha), beta(beta), gamma(gamma) { 
    }
    void optim(optimfn fn, NumericVector init, void * ex) {
      n = init.size();
      nmmin(n, &init[0], &coef[0], &Fmin, fn,
	    &fail, abstol, reltol, ex,
	    alpha, beta, gamma, trace,
	    &fncount, maxit);
    }
    int n, trace, maxit, fail, fncount;
    double abstol, reltol, alpha, beta, gamma, Fmin;
    NumericVector coef;
  };

  // void
  // vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
  //       int maxit, int trace, int *mask,
  //       double abstol, double reltol, int nREPORT, void *ex,
  //       int *fncount, int *grcount, int *fail)

  class BFGS {
  public:
    BFGS(int trace = 0, int maxit = 100, 
	 double abstol = - INFINITY,
	 double reltol = 1.0e-8, int report = 10) : 
      trace(trace), maxit(maxit), report(report), abstol(abstol), reltol(reltol) { }
    void optim(optimfn fn, optimgr gr, NumericVector init, void * ex,
	       double eps = 1.0e-8) {
      n = init.size();
      int mask[n]; for (int i=0; i<n; ++i) mask[i] = 1;
      vmmin(n, &init[0], &Fmin, fn, gr, maxit, trace, mask, abstol, reltol, report,
	    ex, &fncount, &grcount, &fail);
      coef = clone(init);
      hessian = calc_hessian(gr, ex, eps);
    }
    double calc_objective(optimfn fn, NumericVector coef, void * ex) {
      return fn(coef.size(), &coef[0], ex);
    }
    double calc_objective(optimfn fn, void * ex) {
      return fn(coef.size(), &coef[0], ex);
    }
    NumericMatrix calc_hessian(optimgr gr, void * ex, double eps = 1.0e-8) {
      int n = coef.size();
      NumericVector df1(clone(coef));
      NumericVector df2(clone(coef));
      mat hess(n,n,fill::zeros);
      double tmp;
      for(int i=0; i<n; ++i) {
	tmp = coef[i];
	coef[i] += eps;
	gr(n, &coef[0], &df1[0], ex);
	coef[i] = tmp - eps;
	gr(n, &coef[0], &df2[0], ex);
	hess.col(i) = (as<vec>(df1) - as<vec>(df2)) / (2*eps);
	coef[i] = tmp;
      }
      // now symmetrize
      hess = (hess + hess.t()) / 2.0;
      return as<NumericMatrix>(wrap(hess));
    }
    int n, trace, maxit, report, fncount, grcount, fail;
    double abstol, reltol, Fmin;
    NumericVector coef;
    NumericMatrix hessian;
  };

  typedef double (*Brent_fminfn)(double, void *);

  // static
double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		  void *info, double tol)
{
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */

	    d = p / q;
	    u = x + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u, info);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return x;
}

  struct stpm2 {
    mat X, XD, X0;
    vec bhazard,wt,event,wt0;
    int delayed;
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

  template<class Smoother>
  struct pstpm2 {
    NumericVector sp,init;
    mat X, XD, X0;
    vec bhazard,wt,event,wt0;
    int delayed, criterion;
    double reltol;
    std::vector<Smoother> smooth;
  };

  template<class Data>
  double fminfn(int n, double * beta, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vec eta = data->X * vbeta;
    vec h = (data->XD * vbeta) % exp(eta) + data->bhazard;
    double ll = sum(data->wt % data->event % log(h)) - sum(data->wt % exp(eta));
    return -ll;  
  }

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
  }

  template<class Data>
  void grfn(int n, double * beta, double * gr, void *ex) {
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    vec etaD = data->XD * vbeta;
    vec exp_eta = exp(data->X * vbeta);
    // vec h = etaD % exp_eta + data->bhazard;
    vec w23 = data->event % data->wt / (etaD % exp_eta + data->bhazard);
    mat X1 = -data->X;
    X1.each_col() %= (exp_eta % data->wt);
    mat X2 = data->XD;
    X2.each_col() %= (exp_eta % w23);
    mat X3 = data->X;
    X3.each_col() %= (exp_eta % etaD % w23);
    rowvec vgr = -sum(X1+X2+X3,0);
    for (int i = 0; i<n; ++i) gr[i] = vgr[i];
  }

  template<class Smooth>
  void pgrfn(int n, double * beta, double * gr, void *ex) { };

  // Hadamard element-wise multiplication for the _columns_ of a matrix with a vector

  mat rmult(mat m, vec v) {
    mat out(m);
    out.each_col() %= v;
    return out;
  }
  mat lmult(vec v, mat m) {
    return rmult(m,v);
  }

  template<>
  void pgrfn<SmoothLogH>(int n, double * beta, double * gr, void *ex) {
    typedef SmoothLogH Smooth;
    typedef pstpm2<Smooth> Data;
    grfn<Data>(n, beta, gr, ex);
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
    rowvec vgr(n, fill::zeros);
    for (size_t i=0; i < data->smooth.size(); ++i) {
      SmoothLogH smoothi = data->smooth[i];
      vgr.subvec(smoothi.first_para,smoothi.last_para) += 
	(data->sp)[i] * (smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para)).t();
    }
    for (int i = 0; i<n; ++i) gr[i] += vgr[i];
  }

  template<>
  void pgrfn<SmoothHaz>(int n, double * beta, double * gr, void *ex) {
    typedef SmoothHaz Smooth;
    typedef pstpm2<Smooth> Data;
    grfn<Data>(n, beta, gr, ex);
    Data * data = (Data *) ex;
    vec vbeta(beta,n);
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
    for (int i = 0; i<n; ++i) gr[i] += vgr[i];
  }


  RcppExport SEXP optim_stpm2(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
			      SEXP sdelayed, SEXP sX0, SEXP swt0, SEXP sreltol) {

    NumericVector init = as<NumericVector>(sinit);
    int n = init.size();
    NumericMatrix hessian(n,n);

    mat X = as<mat>(sX); 
    mat XD = as<mat>(sXD); 
    vec bhazard = as<vec>(sbhazard);
    vec wt = as<vec>(swt);
    vec event = as<vec>(sevent);

    int delayed = as<int>(sdelayed);
    mat X0(1,1,fill::zeros);
    vec wt0(1,fill::zeros);
    if (delayed == 1) {
      X0 = as<mat>(sX0);
      wt0 = as<vec>(swt0);
    }
    double reltol = as<double>(sreltol);

    stpm2 data = {X, XD, X0, bhazard, wt, event, wt0, delayed};

    BFGS bfgs;
    bfgs.reltol = reltol;
  
    bfgs.optim(fminfn<stpm2>, grfn<stpm2>, init, (void *) &data);

    return List::create(_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));

  }

  template<class Smooth>
  double pstpm2_step_first(double logsp, void * args) {
    typedef pstpm2<Smooth> Data;
    Data * data = (Data *) args;

    data->sp[0] = exp(logsp);

    BFGS bfgs;
    bfgs.reltol = data->reltol;

    bfgs.optim(pfminfn<Smooth>, pgrfn<Smooth>, data->init, (void *) data);

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, (void *) data);

    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));
    double negll = bfgs.calc_objective(fminfn<Data>,(void *) data);
    double gcv =  negll + edf;
    double bic =  negll + edf*log(sum(data->event));
    data->init = bfgs.coef;
    Rprintf("sp=%f\tedf=%f\tnegll=%f\tgcv=%f\tbic=%f\n",data->sp[0],edf,negll,gcv,bic);

    return data->criterion==1 ? gcv : bic;
  }    

  template<class Smooth>
  double pstpm2_step_multivariate(int n, double * logsp, void * args) {
    typedef pstpm2<Smooth> Data;
    Data * data = (Data *) args;

    for (int i=0; i < data->sp.size(); ++i)
      data->sp[i] = bound(exp(logsp[i]),0.001,100.0);

    BFGS bfgs;
    bfgs.reltol = data->reltol;

    bfgs.optim(pfminfn<Smooth>, pgrfn<Smooth>, data->init, (void *) data);

    NumericMatrix hessian0 = bfgs.calc_hessian(grfn<Data>, (void *) data);

    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));
    double negll = bfgs.calc_objective(fminfn<Data>,(void *) data);
    double gcv =  negll + edf;
    double bic =  2.0*negll + edf*log(sum(data->event));
    data->init = bfgs.coef;
    for (int i = 0; i < data->sp.size(); ++i)
      Rprintf("sp[%i]=%f\t",i,data->sp[i]);
    Rprintf("edf=%f\tnegll=%f\tgcv=%f\tbic=%f\n",edf,negll,gcv,bic);

    return data->criterion==1 ? gcv : bic;
  }    

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
  SEXP optim_pstpm2_first(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
			  SEXP sdelayed, SEXP sX0, SEXP swt0,
			  SEXP ssmooth, SEXP ssp, SEXP sreltol_search, SEXP sreltol_final, SEXP salpha, SEXP scriterion) { 

    NumericVector init = as<NumericVector>(sinit);
    NumericVector sp = as<NumericVector>(ssp);
    int n = init.size();
    NumericMatrix hessian(n,n);

    mat X = as<mat>(sX); 
    mat XD = as<mat>(sXD); 
    vec bhazard = as<vec>(sbhazard);
    vec wt = as<vec>(swt);
    vec event = as<vec>(sevent);
    List lsmooth = as<List>(ssmooth);
    int delayed = as<int>(sdelayed);
    double reltol_search = as<double>(sreltol_search);
    double reltol_final = as<double>(sreltol_final);
    double alpha = as<double>(salpha);
    int criterion = as<int>(scriterion);

    mat X0(1,1,fill::zeros);
    vec wt0(1,fill::zeros);
    if (delayed == 1) {
      X0 = as<mat>(sX0);
      wt0 = as<vec>(swt0);
    }

    std::vector<Smooth> smooth = read_smoothers<Smooth>(lsmooth);

    pstpm2<Smooth> data = {sp, init, X, XD, X0, bhazard, wt, event, wt0, delayed, criterion, reltol_search, smooth};

    double opt_sp = alpha * exp(Brent_fmin(log(0.001),log(100.0),&pstpm2_step_first<Smooth>,&data,1.0e-2));
    data.sp[0] = opt_sp;

    BFGS bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = reltol_final;
    bfgs.optim(pfminfn<Smooth>, pgrfn<Smooth>, data.init, (void *) &data);

    hessian = bfgs.calc_hessian(pgrfn<Smooth>, (void *) &data);

    return List::create(_("sp")=wrap(opt_sp),
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(hessian));

  }

  RcppExport SEXP optim_pstpm2LogH_first(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, 
					 SEXP swt, SEXP sevent,
					 SEXP sdelayed, SEXP sX0, SEXP swt0,
					 SEXP ssmooth, SEXP ssp, SEXP sreltol_search, 
					 SEXP sreltol_final, SEXP salpha, SEXP scriterion) {
    return optim_pstpm2_first<SmoothLogH>(sinit, sX, sXD, sbhazard, swt, sevent,
					  sdelayed, sX0, swt0,
					  ssmooth, ssp, sreltol_search, sreltol_final, salpha, scriterion);
  }

  RcppExport SEXP optim_pstpm2Haz_first(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, 
					 SEXP swt, SEXP sevent,
					 SEXP sdelayed, SEXP sX0, SEXP swt0,
					 SEXP ssmooth, SEXP ssp, SEXP sreltol_search, 
					SEXP sreltol_final, SEXP salpha, SEXP scriterion) {
    return optim_pstpm2_first<SmoothHaz>(sinit, sX, sXD, sbhazard, swt, sevent,
					  sdelayed, sX0, swt0,
					 ssmooth, ssp, sreltol_search, sreltol_final, salpha, scriterion);
  }

  template<class Smooth>
  SEXP optim_pstpm2_multivariate(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
				 SEXP sdelayed, SEXP sX0, SEXP swt0,
				 SEXP ssmooth, SEXP ssp, SEXP sreltol_search, SEXP sreltol_final, SEXP salpha, SEXP scriterion) {

    NumericVector init = as<NumericVector>(sinit);
    NumericVector sp = as<NumericVector>(ssp);
    // NumericVector x(init.size()); 
    int n = init.size();
    NumericMatrix hessian(n,n);

    mat X = as<mat>(sX); 
    mat XD = as<mat>(sXD); 
    vec bhazard = as<vec>(sbhazard);
    vec wt = as<vec>(swt);
    vec event = as<vec>(sevent);
    List lsmooth = as<List>(ssmooth);
    int delayed = as<int>(sdelayed);
    int criterion = as<int>(scriterion);
    double reltol_search = as<double>(sreltol_search);
    double reltol_final = as<double>(sreltol_final);
    double alpha = as<double>(salpha);

    mat X0(1,1,fill::zeros);
    vec wt0(1,fill::zeros);
    if (delayed == 1) {
      X0 = as<mat>(sX0);
      wt0 = as<vec>(swt0);
    }

    std::vector<Smooth> smooth = read_smoothers<Smooth>(lsmooth);

    pstpm2<Smooth> data = {sp, init, X, XD, X0, bhazard, wt, event, wt0, delayed, criterion, reltol_search, smooth};

    NelderMead nm;
    nm.reltol = 1.0e-4;
    nm.maxit = 500;
    nm.optim(pstpm2_step_multivariate<Smooth>, data.sp, (void *) &data);

    for (int i=0; i < nm.coef.size(); ++i)
      data.sp[i] = alpha * exp(nm.coef[i]); // same alpha for all smoothers??

    BFGS bfgs;
    bfgs.coef = data.init;
    bfgs.reltol = reltol_final;
    bfgs.optim(pfminfn<Smooth>, pgrfn<Smooth>, data.init, (void *) &data);

    hessian = bfgs.calc_hessian(pgrfn<Smooth>, (void *) &data);

    return List::create(_("sp")=wrap(data.sp),
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(hessian));

  }

  RcppExport SEXP optim_pstpm2LogH_multivariate(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
				 SEXP sdelayed, SEXP sX0, SEXP swt0,
						SEXP ssmooth, SEXP ssp, SEXP sreltol_search, SEXP sreltol_final, SEXP salpha, SEXP scriterion) {
    return optim_pstpm2_multivariate<SmoothLogH>(sinit, sX, sXD, sbhazard, swt, sevent,
				 sdelayed, sX0, swt0,
						 ssmooth, ssp, sreltol_search, sreltol_final, salpha, scriterion);
  }

  RcppExport SEXP optim_pstpm2Haz_multivariate(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
				 SEXP sdelayed, SEXP sX0, SEXP swt0,
					       SEXP ssmooth, SEXP ssp, SEXP sreltol_search, SEXP sreltol_final, SEXP salpha, SEXP scriterion) {
    return optim_pstpm2_multivariate<SmoothHaz>(sinit, sX, sXD, sbhazard, swt, sevent,
				 sdelayed, sX0, swt0,
						ssmooth, ssp, sreltol_search, sreltol_final, salpha, scriterion);
  }


  template<class Smooth>
  SEXP optim_pstpm2_fixedsp(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
  		  SEXP sdelayed, SEXP sX0, SEXP swt0,
			  SEXP ssmooth, SEXP ssp, SEXP sreltol) { 

    NumericVector init = as<NumericVector>(sinit);
    NumericVector sp = as<NumericVector>(ssp);
    int n = init.size();
    NumericMatrix hessian(n,n);

    mat X = as<mat>(sX); 
    mat XD = as<mat>(sXD); 
    vec bhazard = as<vec>(sbhazard);
    vec wt = as<vec>(swt);
    vec event = as<vec>(sevent);
    List lsmooth = as<List>(ssmooth);
    int delayed = as<int>(sdelayed);
    double reltol = as<double>(sreltol);

    mat X0(1,1,fill::zeros);
    vec wt0(1,fill::zeros);
    if (delayed == 1) {
      X0 = as<mat>(sX0);
      wt0 = as<vec>(swt0);
    }

    std::vector<Smooth> smooth = read_smoothers<Smooth>(lsmooth);

    pstpm2<Smooth> data = {sp, init, X, XD, X0, bhazard, wt, event, wt0, delayed, 1, reltol, smooth};

    BFGS bfgs;
    bfgs.coef = init;
    bfgs.reltol = reltol;
    bfgs.optim(pfminfn<Smooth>, pgrfn<Smooth>, data.init, (void *) &data);

    hessian = bfgs.calc_hessian(pgrfn<Smooth>, (void *) &data);

    return List::create(_("sp")=wrap(sp),
			_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(hessian));

  }

  RcppExport SEXP optim_pstpm2LogH_fixedsp(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, 
					 SEXP swt, SEXP sevent,
					 SEXP sdelayed, SEXP sX0, SEXP swt0,
					 SEXP ssmooth, SEXP ssp, 
					 SEXP sreltol) {
    return optim_pstpm2_fixedsp<SmoothLogH>(sinit, sX, sXD, sbhazard, swt, sevent,
					  sdelayed, sX0, swt0,
					  ssmooth, ssp, sreltol);
  }

  RcppExport SEXP optim_pstpm2Haz_fixedsp(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, 
					 SEXP swt, SEXP sevent,
					 SEXP sdelayed, SEXP sX0, SEXP swt0,
					 SEXP ssmooth, SEXP ssp, 
					 SEXP sreltol) {
    return optim_pstpm2_fixedsp<SmoothHaz>(sinit, sX, sXD, sbhazard, swt, sevent,
					  sdelayed, sX0, swt0,
					  ssmooth, ssp, sreltol);
  }



  template<class Smooth>
  SEXP test_pstpm2(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
  		  SEXP sdelayed, SEXP sX0, SEXP swt0,
			  SEXP ssmooth, SEXP ssp, SEXP sreltol) { 

    NumericVector init = as<NumericVector>(sinit);
    NumericVector sp = as<NumericVector>(ssp);
    int n = init.size();
    NumericMatrix hessian(n,n);

    mat X = as<mat>(sX); 
    mat XD = as<mat>(sXD); 
    vec bhazard = as<vec>(sbhazard);
    vec wt = as<vec>(swt);
    vec event = as<vec>(sevent);
    List lsmooth = as<List>(ssmooth);
    int delayed = as<int>(sdelayed);
    double reltol = as<double>(sreltol);

    mat X0(1,1,fill::zeros);
    vec wt0(1,fill::zeros);
    if (delayed == 1) {
      X0 = as<mat>(sX0);
      wt0 = as<vec>(swt0);
    }

    std::vector<Smooth> smooth = read_smoothers<Smooth>(lsmooth);

    pstpm2<Smooth> data = {sp, init, X, XD, X0, bhazard, wt, event, wt0, delayed, 1, reltol, smooth};

    double pfmin = pfminfn<Smooth>(n,&init[0],(void *) &data);
    NumericVector gr(n);
    pgrfn<Smooth>(n,&init[0], &gr[0], (void *) &data);

    return List::create(_("sp")=wrap(sp),
			_("pfmin")=wrap(pfmin),
			_("gr")=wrap(gr)
			);

  }

  RcppExport SEXP test_pstpm2LogH(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, 
					 SEXP swt, SEXP sevent,
					 SEXP sdelayed, SEXP sX0, SEXP swt0,
					 SEXP ssmooth, SEXP ssp, 
					 SEXP sreltol) {
    return test_pstpm2<SmoothLogH>(sinit, sX, sXD, sbhazard, swt, sevent,
					  sdelayed, sX0, swt0,
					  ssmooth, ssp, sreltol);
  }

  RcppExport SEXP test_pstpm2Haz(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, 
					 SEXP swt, SEXP sevent,
					 SEXP sdelayed, SEXP sX0, SEXP swt0,
					 SEXP ssmooth, SEXP ssp, 
					 SEXP sreltol) {
    return test_pstpm2<SmoothHaz>(sinit, sX, sXD, sbhazard, swt, sevent,
					  sdelayed, sX0, swt0,
					  ssmooth, ssp, sreltol);
  }



  // R CMD INSTALL ~/src/R/microsimulation
  // R -q -e "require(microsimulation); .Call('test_nmmin',1:2,PACKAGE='microsimulation')"

  // .Call("optim_stpm2",init,X,XD,rep(bhazard,nrow(X)),wt,ifelse(event,1,0),package="rstpm2")

} // anonymous namespace
