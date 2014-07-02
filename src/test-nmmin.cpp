#include <R_ext/Applic.h>
#include <RcppArmadillo.h>
#include <vector>
#include <float.h> /* DBL_EPSILON */

namespace {

  using namespace Rcpp;
  using namespace arma;

  typedef double optimfn(int, double *, void *);
  typedef void optimgr(int, double *, double *, void *);

  void Rprint(NumericMatrix m) {
    for (int i=0; i<m.nrow(); ++i) {
      for (int j=0; j<m.ncol(); ++j) 
	Rprintf("%f ", m(i,j));
      Rprintf("\n");
    }
  }

  void Rprint(mat m) {
    for (int i=0; i<m.n_rows; ++i) {
      for (int j=0; j<m.n_cols; ++j) 
	Rprintf("%f ", m(i,j));
      Rprintf("\n");
    }
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
    void optim(optimfn fn, NumericVector init, NumericVector x, void * ex) {
      n = init.size();
      nmmin(n, &init[0], &x[0], &Fmin, fn,
	    &fail, abstol, reltol, ex,
	    alpha, beta, gamma, trace,
	    &fncount, maxit);
      coef = clone(init);
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
      trace(trace), maxit(maxit), abstol(abstol), reltol(reltol), report(report) { }
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

  struct Smooth {
    int first_para, last_para;
    mat S;
  };

  struct pstpm2 {
    NumericVector sp,init;
    mat X, XD, X0;
    vec bhazard,wt,event,wt0;
    int delayed;
    std::vector<Smooth> smooth;
  };

  double fminfn(int n, double * beta, void *ex) {
    stpm2 * data = (stpm2 *) ex;
    vec vbeta(beta,n);
    vec eta = data->X * vbeta;
    vec h = (data->XD * vbeta) % exp(eta) + data->bhazard;
    double ll = sum(data->wt % data->event % log(h)) - sum(data->wt % exp(eta));
    // for (int i = 0; i < n; ++i) 
    //   Rprintf("%f ",vbeta[i]);
    // Rprintf("\nll=% f\n",ll);
    return -ll;  
  }

  void grfn(int n, double * beta, double * gr, void *ex) {
    stpm2 * data = (stpm2 *) ex;
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

  double pfminfn(int n, double * beta, void *ex) {
    pstpm2 * data = (pstpm2 *) ex;
    vec vbeta(beta,n);
    vec eta = data->X * vbeta;
    vec h = (data->XD * vbeta) % exp(eta) + data->bhazard;
    double ll = sum(data->wt % data->event % log(h)) - sum(data->wt % exp(eta));
    for (int i=0; i < data->smooth.size(); ++i) {
      Smooth smoothi = data->smooth[i];
      ll -= (data->sp)[i]/2 * 
	dot(vbeta.subvec(smoothi.first_para,smoothi.last_para),
	    smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para));
    }
    return -ll;  
  }


  void pgrfn(int n, double * beta, double * gr, void *ex) {
    pstpm2 * data = (pstpm2 *) ex;
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
    for (int i=0; i < data->smooth.size(); ++i) {
      Smooth smoothi = data->smooth[i];
      vgr.subvec(smoothi.first_para,smoothi.last_para) += 
	(data->sp)[i] * (smoothi.S * vbeta.subvec(smoothi.first_para,smoothi.last_para)).t();
    }
    for (int i = 0; i<n; ++i) gr[i] = vgr[i];
  }

  void pgrfn_unpenalised(int n, double * beta, double * gr, void *ex) {
    pstpm2 * data = (pstpm2 *) ex;
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


  RcppExport SEXP optim_stpm2(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
			      SEXP sdelayed, SEXP sX0, SEXP swt0) {

    NumericVector init = as<NumericVector>(sinit);
    // NumericVector x(init.size()); 
    NumericVector x2(clone(init));
    int n = x2.size();
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

    stpm2 data = {X, XD, X0, bhazard, wt, event, wt0, delayed};

    // NelderMead nm;
    // nm.maxit = 1000;
    // nm.reltol = 1.0e-10;
    BFGS bfgs;
    bfgs.reltol = 1.0e-10;
  
    // nm.optim(fminfn, init, x, (void *) &data);

    // Rprintf("beta[0]=%g\nbeta[1]=%g\n",x[0],x[1]);
    // Rprintf("fail=%i\nFmin=%g\nfncount=%i\n",nm.fail,nm.Fmin,nm.fncount);

    bfgs.optim(fminfn, grfn, init, (void *) &data);

    // Rprintf("beta[0]=%g\nbeta[1]=%g\n",x2[0],x2[1]);
    // Rprintf("fail=%i\nFmin=%g\nfncount=%i\ngrcount=%i\n",bfgs.fail,bfgs.Fmin,bfgs.fncount,bfgs.grcount);

    return List::create(_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian));

  }

  double pstpm2_step(double logsp, void * args) {
    pstpm2 * data = (pstpm2 *) args;

    data->sp[0] = exp(logsp);

    BFGS bfgs;
    bfgs.reltol = 1.0e-10;

    bfgs.optim(pfminfn, pgrfn, data->init, (void *) data);

    NumericMatrix hessian0 = bfgs.calc_hessian(pgrfn_unpenalised, (void *) data);

    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));
    double negll = bfgs.calc_objective(pfminfn,(void *) data);
    double gcv =  negll + edf;
    data->init = bfgs.coef;
    Rprintf("sp=%f\tedf=%f\tnegll=%f\tgcv=%f\n",data->sp[0],edf,negll,gcv);

    return gcv;
  }    

  RcppExport SEXP optim_pstpm2(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
			       SEXP sdelayed, SEXP sX0, SEXP swt0,
			       SEXP ssmooth, SEXP ssp) {

    NumericVector init = as<NumericVector>(sinit);
    NumericVector sp = as<NumericVector>(ssp);
    // NumericVector x(init.size()); 
    NumericVector x2(clone(init));
    int n = init.size();
    NumericMatrix hessian(n,n);

    mat X = as<mat>(sX); 
    mat XD = as<mat>(sXD); 
    vec bhazard = as<vec>(sbhazard);
    vec wt = as<vec>(swt);
    vec event = as<vec>(sevent);
    List lsmooth = as<List>(ssmooth);
    int delayed = as<int>(sdelayed);

    mat X0(1,1,fill::zeros);
    vec wt0(1,fill::zeros);
    if (delayed == 1) {
      X0 = as<mat>(sX0);
      wt0 = as<vec>(swt0);
    }

    std::vector<Smooth> smooth;
    for(int i=0; i<lsmooth.size(); ++i) {
      List lsmoothi = as<List>(lsmooth[i]);
      List lsmoothiS = as<List>(lsmoothi["S"]);
      Smooth smoothi =  {
	as<int>(lsmoothi["first.para"]) - 1, 
	as<int>(lsmoothi["last.para"]) - 1, 
	as<mat>(lsmoothiS[0])
      };
      smooth.push_back(smoothi);
    }
    pstpm2 data = {sp, init, X, XD, X0, bhazard, wt, event, wt0, delayed, smooth};

    BFGS bfgs;
    bfgs.reltol = 1.0e-10;

    // Rprintf("Objective(init)=%f.\n",bfgs.calc_objective(pfminfn,data.init, (void *) &data));

    bfgs.optim(pfminfn, pgrfn, data.init, (void *) &data);

    NumericMatrix hessian0 = bfgs.calc_hessian(pgrfn_unpenalised, (void *) &data);

    double edf = trace(solve(as<mat>(bfgs.hessian),as<mat>(hessian0)));
    double gcv = bfgs.Fmin - edf;

    return List::create(_("coef")=wrap(bfgs.coef),
			_("hessian")=wrap(bfgs.hessian),
			_("hessian0")=wrap(hessian0),
			_("ll")=wrap(-bfgs.Fmin),
			_("gcv") = wrap(gcv));

  }

  RcppExport SEXP optim_pstpm2_gcv(SEXP sinit, SEXP sX, SEXP sXD, SEXP sbhazard, SEXP swt, SEXP sevent,
			       SEXP sdelayed, SEXP sX0, SEXP swt0,
			       SEXP ssmooth, SEXP ssp) {

    NumericVector init = as<NumericVector>(sinit);
    NumericVector sp = as<NumericVector>(ssp);
    // NumericVector x(init.size()); 
    NumericVector x2(clone(init));
    int n = init.size();
    NumericMatrix hessian(n,n);

    mat X = as<mat>(sX); 
    mat XD = as<mat>(sXD); 
    vec bhazard = as<vec>(sbhazard);
    vec wt = as<vec>(swt);
    vec event = as<vec>(sevent);
    List lsmooth = as<List>(ssmooth);
    int delayed = as<int>(sdelayed);

    mat X0(1,1,fill::zeros);
    vec wt0(1,fill::zeros);
    if (delayed == 1) {
      X0 = as<mat>(sX0);
      wt0 = as<vec>(swt0);
    }

    std::vector<Smooth> smooth;
    for(int i=0; i<lsmooth.size(); ++i) {
      List lsmoothi = as<List>(lsmooth[i]);
      List lsmoothiS = as<List>(lsmoothi["S"]);
      Smooth smoothi =  {
	as<int>(lsmoothi["first.para"]) - 1, 
	as<int>(lsmoothi["last.para"]) - 1, 
	as<mat>(lsmoothiS[0])
      };
      smooth.push_back(smoothi);
    }

    pstpm2 data = {sp, init, X, XD, X0, bhazard, wt, event, wt0, delayed, smooth};

    double opt_sp = exp(Brent_fmin(log(0.001),log(10.0),&pstpm2_step,&data,1.0e-2));

    BFGS bfgs;
    bfgs.coef = data.init;

    hessian = bfgs.calc_hessian(pgrfn, (void *) &data);

    return List::create(_("sp")=wrap(opt_sp),
			_("coef")=wrap(data.init),
			_("hessian")=wrap(hessian));

  }


  // R CMD INSTALL ~/src/R/microsimulation
  // R -q -e "require(microsimulation); .Call('test_nmmin',1:2,PACKAGE='microsimulation')"

  // .Call("optim_stpm2",init,X,XD,rep(bhazard,nrow(X)),wt,ifelse(event,1,0),package="rstpm2")

} // anonymous namespace
