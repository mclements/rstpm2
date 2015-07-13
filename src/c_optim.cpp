#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <vector>
#include <map>
#include <float.h> /* DBL_EPSILON */
#include "c_optim.h"

// #include "uncmin.cpp"

namespace rstpm2 {

  using namespace Rcpp;

  double min(double a, double b) { return a < b ? a : b; }
  double max(double a, double b) { return a < b ? b : a; }
  double bound(double x, double lower, double upper) { return x < lower ? lower : (x > upper ? upper : x); }

  // void nmmin(int n, double *Bvec, double *X, double *Fmin, optimfn fminfn,
  // 	   int *fail, double abstol, double intol, void *ex,
  // 	   double alpha, double bet, double gamm, int trace,
  // 	   int *fncount, int maxit)

  NelderMead::NelderMead(int trace, int maxit, 
			 double abstol, double reltol, 
			 double alpha, double beta, double gamma,
			 double epshess, bool hessianp) : 
      trace(trace), maxit(maxit), abstol(abstol), reltol(reltol), 
      alpha(alpha), beta(beta), gamma(gamma), epshess(epshess), hessianp(hessianp) { 
    }
  void NelderMead::optim(optimfn fn, NumericVector init, void * ex) {
    n = init.size();
    coef = clone(init);
    nmmin(n, &init[0], &coef[0], &Fmin, fn,
	    &fail, abstol, reltol, ex,
	    alpha, beta, gamma, trace,
	    &fncount, maxit);
      if (hessianp)
	hessian = calc_hessian(fn, ex);
    }
  NumericMatrix NelderMead::calc_hessian(optimfn fn, void * ex) {
      int n = coef.size();
      NumericMatrix hess(n,n);
      double tmpi,tmpj,f1,f0,fm1,hi,hj,fij,fimj,fmij,fmimj;
      f0 = fn(n,&coef[0],ex);
      for(int i=0; i<n; ++i) {
	tmpi = coef[i];
	hi = epshess*(1.0+abs(tmpi));
	coef[i] = tmpi + hi;
	f1=fn(n, &coef[0], ex);
	coef[i] = tmpi - hi;
	fm1=fn(n, &coef[0], ex);
	// hess(i,i) = (-f2 +16.0*f1 - 30.0*f0 + 16.0*fm1 -fm2)/(12.0*hi*hi);
	hess(i,i) = (f1 - 2.0*f0 + fm1)/(hi*hi);
	coef[i] = tmpi;
	for (int j=i; j<n; ++j) {
	  if (i != j) {
	    tmpj = coef[j];
	    hj = epshess*(1.0+abs(tmpj));
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
	    hess(j,i) = hess(i,j) = (fij-fimj-fmij+fmimj)/(4.0*hi*hj);
	    coef[i] = tmpi;
	    coef[j] = tmpj;
	  }
	}
      }
      return wrap(hess);
    }

  // void
  // vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
  //       int maxit, int trace, int *mask,
  //       double abstol, double reltol, int nREPORT, void *ex,
  //       int *fncount, int *grcount, int *fail)
  BFGS::BFGS(int trace, int maxit, 
	 double abstol,
	     double reltol, int report, double epshess, bool hessianp) : 
    trace(trace), maxit(maxit), report(report), abstol(abstol), reltol(reltol), epshess(epshess), hessianp(hessianp) { }
  void BFGS::optim(optimfn fn, optimgr gr, NumericVector init, void * ex) {
      n = init.size();
      std::vector<int> mask(n,1); 
      vmmin(n, &init[0], &Fmin, fn, gr, maxit, trace, &mask[0], abstol, reltol, report,
	    ex, &fncount, &grcount, &fail);
      coef = clone(init);
      if (hessianp)
	hessian = calc_hessian(gr, ex);
    }
  double BFGS::calc_objective(optimfn fn, NumericVector coef, void * ex) {
      return fn(coef.size(), &coef[0], ex);
    }
  double BFGS::calc_objective(optimfn fn, void * ex) {
      return fn(coef.size(), &coef[0], ex);
    }
  NumericMatrix BFGS::calc_hessian(optimgr gr, void * ex) {
      int n = coef.size();
      NumericVector df1(clone(coef));
      NumericVector df2(clone(coef));
      NumericMatrix hess(n,n);
      double tmp;
      for(int i=0; i<n; ++i) {
	tmp = coef[i];
	coef[i] += epshess;
	gr(n, &coef[0], &df1[0], ex);
	coef[i] = tmp - epshess;
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
      return wrap(hess);
    }


// void
// optif9(int nr, int n, double *x, fcn_p fcn, fcn_p d1fcn, d2fcn_p d2fcn,
//        void *state, double *typsiz, double fscale, int method,
//        int iexp, int *msg, int ndigit, int itnlim, int iagflg, int iahflg,
//        double dlt, double gradtl, double stepmx, double steptl,
//        double *xpls, double *fpls, double *gpls, int *itrmcd, double *a,
//        double *wrk, int *itncnt)
  Nlm::Nlm(double fscale,
	   int method,
	   int iexp,
	   int msg,
	   int ndigit,
	   int itnlim,
	   int iagflg,
	   int iahflg,
	   double dlt,
	   double gradtl, // cf. epshess
	   double stepmx,
	   double steptl,
	   int itrmcd,
	   int itncnt,
	   bool hessianp
	   ) : fscale(fscale), method(method), iexp(iexp), msg(msg),
	       ndigit(ndigit), itnlim(itnlim), iagflg(iagflg), 
	       iahflg(iahflg), dlt(dlt), gradtl(gradtl), stepmx(stepmx),
	       steptl(steptl), itrmcd(itrmcd), itncnt(itncnt), hessianp(hessianp) { }
  void Nlm::optim(fcn_p fcn, fcn_p d1fcn, NumericVector init, void * state) {
      int n;
      n = init.size();
      std::vector<double> typsize(n,1.0), gpls(n,0.0), a(n*n,0.0), wrk(n*8,0.0);
      double norm, fpls;
      NumericVector xpls(n);
      // stepmax calculations
      if (stepmx == -1.0) {
	norm = 0.0;
	for (int i=0; i<n; ++i)
	  norm += init[i]*init[i]/typsize[i]/typsize[i];
	norm = sqrt(norm);
	stepmx = norm < 1.0 ? 1000.0 : norm*1000.0;
      }
      // call the optimizer
      optif9(n, n, &init[0], fcn, d1fcn, (d2fcn_p) 0, state, &typsize[0], fscale, method, 
	     iexp, &msg, ndigit, itnlim, iagflg, iahflg,
	     dlt, gradtl, stepmx, steptl,
	     &xpls[0], &fpls, &gpls[0], &itrmcd, &a[0],
	     &wrk[0], &itncnt);
      coef = clone(xpls);
      if (hessianp)
	hessian = calc_hessian(d1fcn, state);
    }
  double Nlm::calc_objective(fcn_p fn, NumericVector coef, void * ex) {
    double f;
    fn(coef.size(), &coef[0], &f, ex);
    return f;
  }
  double Nlm::calc_objective(fcn_p fn, void * ex) {
    double f;
    fn(coef.size(), &coef[0], &f, ex);
    return f;
  }
  NumericMatrix Nlm::calc_hessian(fcn_p gr, void * ex) {
      int n = coef.size();
      NumericVector df1(clone(coef));
      NumericVector df2(clone(coef));
      NumericMatrix hess(n,n);
      double tmp;
      for(int i=0; i<n; ++i) {
	tmp = coef[i];
	coef[i] += gradtl;
	gr(n, &coef[0], &df1[0], ex);
	coef[i] = tmp - gradtl;
	gr(n, &coef[0], &df2[0], ex);
	for (int j=i; j<n; ++j)
	  hess(j,i) = hess(i,j) = (df1[j] - df2[j]) / (2*gradtl);
	coef[i] = tmp;
      }
      return wrap(hess);
  }
  void Nlm::set_print_level(int print_level) {
    if (print_level == 0) msg = 9;
    if (print_level == 1) msg = 1;
    if (print_level >= 2) msg = 17;
  }

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

  // R CMD INSTALL ~/src/R/microsimulation
  // R -q -e "require(microsimulation); .Call('test_nmmin',1:2,PACKAGE='microsimulation')"

  // .Call("optim_stpm2",init,X,XD,rep(bhazard,nrow(X)),wt,ifelse(event,1,0),package="rstpm2")

} // namespace rstpm2
