#include <math.h>
#include <float.h>
#include <Rmath.h> /* for fmax2, fmin2, imin2 */
#include <RcppArmadillo.h>
#include <vector>

/*
  - Adapted from src/appl/integrate.c, which is a modified f2c translation of QUADPACK.
  - Principle: change as little of the code as possible
  - Re-factor to carefully use vec and mat -- but use C arrays if possible and loop when required
  - f2c uses 1-based indexes -- the re-factor uses some 0-based indexes (elist, rlist, elistSum)
  - As a reminder, function f: R^ny -> R^ny for each integration point and ny function inputs/outputs
    (whereas integrate.c has f: R^n -> R^n for n integration points)
  - This now uses templates, which allows for vectorised integration in C++ as well as in R
  - Potential issue: use of any() or all() in the conditions may be incorrect. 
    Principle: we want to be conservative.
  - Potential issue: limit*ny and 52*ny may be large (leading ot C stack overflow):(
  - Limitation: in some cases, this is not much faster than using integrate() in a for loop in R
  - TODO: can this be extended to allow for vectorised integration in C++?

  - GPL>=3
  - Mark Clements 2023-07-31
*/

//[[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using std::vector;

template<typename F>
vec eval_f(const F f, const vec x) {
  return f(x);
}
template<>
vec eval_f(const Rcpp::Function f, const vec x) {
  return Rcpp::as<vec>(Rcpp::wrap(f(Rcpp::wrap(x))));
}

template<typename F>
void vrdqk21(const F f,
	     const vec, const vec,
	     const double, const double, vec *, vec *, vec *, vec *);

template<typename F>
void vrdqk15i(const F f,
	      const vec boun,
	      const int inf, double *a, double *b,
	      vec *result,
	      vec *abserr, vec *resabs, vec *resasc);

static void rdqpsrt(const int, int *, int *, double *, double *, int *, int *);

static void rdqelg(int *, double *, double *, double *, double *, int *);

template<typename F>
void vrdqagse(const F f, const vec a, const vec b, const double 
	      epsabs, const double epsrel, const int limit, const int ny,
	      double *result, double *abserr, int *neval, int *ier, double *alist,
	      double *blist, double *rlist, double *elist, int *
	      iord, int *last);

template<typename F>
void vrdqagie(const F f, const vec bound, const int inf,
	      const double epsabs, const double epsrel, const int limit, const int ny,
	      double *resultp, double *abserrp, int *neval, int *ier, double *alist,
	      double *blist, double *rlistp, double *elistp, int *iord, int *last);

template<typename F>
void vRdqagi(const F f, const vec bound, const int inf,
	     const double epsabs, const double epsrel, const int limit, const int ny,
	     double *result, double *abserr, int *neval, int *ier,
	     int *lenw, int *last,
	     int *iwork, double *work)
{
  int l1, l2, l3;
  *ier = 6;
  *neval = 0;
  *last = 0;
  for (int i=0; i<ny; i++) {
    result[i] = 0.0;
    abserr[i] = 0.0;
  }
  if (limit < 1 || *lenw < limit*2 + limit*ny*2) return;

  /*         prepare call for vrdqagie. */

  // double work[2*limit + 2*limit*ny];
  l1 = limit;
  l2 = limit + l1;
  l3 = limit*ny + l2;

  vrdqagie<F>(f, bound, inf, epsabs, epsrel, limit, ny,
	   result, abserr, neval, ier,
	   work, &work[l1], &work[l2], &work[l3], iwork, last);

  return;
} /* vRdqags_ */

template<typename F>
void vRdqags(const F f, const vec a, const vec b,
	     const double epsabs, const double epsrel, const int ny,
	     double *result, double *abserr, int *neval, int *ier,
	     const int limit, int *lenw, int *last, int *iwork, double *work)
{
  int l1, l2, l3;

  /*         check validity of limit and lenw. */

  *ier = 6;
  *neval = 0;
  *last = 0;
  for (int i=0; i<ny; i++) {
    result[i] = 0.0;
    abserr[i] = 0.0;
  }
  if (limit < 1 || *lenw < limit*2 + limit*ny*2) return;

  /*         prepare call for dqagse. */

  // double work[2*limit + 2*limit*ny];
  l1 = limit;
  l2 = limit + l1;
  l3 = limit*ny + l2;

  vrdqagse<F>(f, a, b, epsabs, epsrel, limit, ny, result, abserr, neval, ier,
	      work, &work[l1], &work[l2], &work[l3], iwork, last);

  return;
} /* vRdqags_ */

template<typename F>
void vrdqagse(const F f, const vec lower, const vec upper, const double 
	      epsabs, const double epsrel, const int limit, const int ny,
	      double *resultp, double *abserrp, int *neval, int *ier, double *alist,
	      double *blist, double *rlistp, double *elistp, int *
	      iord, int *last)
{
  /* Local variables */
  Rboolean noext, extrap;
  int k,ksgn, nres;
  int ierro;
  int ktmin, nrmax;
  int iroff1, iroff2, iroff3;
  int id;
  int numrl2;
  int jupbnd;
  int maxerr;
  double *res3la = R_Calloc(3*ny, double);;
  double *rlist2 = R_Calloc(52*ny, double);
  vec abseps= zeros(ny), area1= zeros(ny), area2= zeros(ny), area12= zeros(ny);
  double epmach;
  double a, b, a1, a2, b1, b2, oflow, uflow;
  vec defab1= zeros(ny), defab2= zeros(ny), reseps= zeros(ny);
  vec area= zeros(ny), defabs= zeros(ny), resabs= zeros(ny), dres= zeros(ny), errbnd= zeros(ny);
  vec error1= zeros(ny), error2= zeros(ny), erro12= zeros(ny), correc= zeros(ny), erlarg= zeros(ny), errsum= zeros(ny);
  vec errmax= zeros(ny), erlast= zeros(ny), ertest= zeros(ny);
  double errmaxsum, small = 0.0;

  // vec result(resultp, limit, false); 
  // vec abserr(abserrp, limit, false); 
  mat elist(elistp, ny, limit, false);
  mat rlist(rlistp, ny, limit, false);
  double *elistSum = R_Calloc(limit, double);;
    
  /* Parameter adjustments */
  --iord;
  // --elist;
  // --rlist;
  --blist;
  --alist;

  /* Function Body */
  epmach = DBL_EPSILON;

  /*            test on validity of parameters */
  /*            ------------------------------ */
  *ier = 0;
  *neval = 0;
  *last = 0;
  vec result = zeros(ny);
  vec abserr = zeros(ny);
  a = 0.0;
  b = 1.0;
  alist[1] = a;
  blist[1] = b;
  rlist.col(0) = zeros(ny);
  elist.col(0) = zeros(ny);
  elistSum[0] = 0.0;
  if (epsabs <= 0. && epsrel < R::fmax2(epmach * 50., 5e-29)) {
    *ier = 6;
    return;
  }

  /*           first approximation to the integral */
  /*           ----------------------------------- */

  uflow = DBL_MIN;
  oflow = DBL_MAX;
  ierro = 0;
  vrdqk21<F>(f, lower, upper, a, b, &result, &abserr, &defabs, &resabs);

  /*           test on accuracy. */

  dres = abs(result);
  errbnd = max(ones(ny)*epsabs, epsrel * dres);
  *last = 1;
  rlist.col(0) = result;
  elist.col(0) = abserr;
  elistSum[0] = max(abserr);
  iord[1] = 1;
  if (all(abserr <= epmach * 100. * defabs && abserr > errbnd)) // "roundoff error was detected"
    *ier = 2;
  if (limit == 1) // "maximum number of subdivisions reached"
    *ier = 1;
  if (*ier != 0
      || all(abserr <= errbnd && abserr != resabs) // "OK"
      || all(abserr == 0.)) // "OK"
    goto L140;

  /*           initialization */
  /*           -------------- */

  for (int i=0; i<ny; ++i)
    rlist2[i*52] = result[i];
  errmax = abserr;
  maxerr = 1; // base 1
  area = result;
  errsum = abserr;
  for (int i=0; i<ny; ++i)
    abserr[i] = oflow;
  nrmax = 1;
  nres = 0;
  numrl2 = 2;
  ktmin = 0;
  extrap = FALSE;
  noext = FALSE;
  iroff1 = 0;
  iroff2 = 0;
  iroff3 = 0;
  ksgn = -1;
  if (any(dres >= (1. - epmach * 50.) * defabs)) {
    ksgn = 1;
  }

  /*           main do-loop */
  /*           ------------ */

  for (*last = 2; *last <= limit; ++(*last)) {

    /*           bisect the subinterval with the nrmax-th largest error estimate. */

    a1 = alist[maxerr];
    b1 = (alist[maxerr] + blist[maxerr]) * .5;
    a2 = b1;
    b2 = blist[maxerr];
    erlast = errmax;
    vrdqk21<F>(f, lower, upper, a1, b1, &area1, &error1, &resabs, &defab1);
    vrdqk21<F>(f, lower, upper, a2, b2, &area2, &error2, &resabs, &defab2);

    /*           improve previous approximations to integral
		 and error and test for accuracy. */

    area12 = area1 + area2;
    erro12 = error1 + error2;
    errsum = errsum + erro12 - errmax;
    area = area + area12 - rlist.col(maxerr-1);
    if (all(defab1 != error1 && defab2 != error2)) { // *any* or all?

      if (all(abs(rlist.col(maxerr-1) - area12) <= abs(area12) * 1e-5 && // any or all?
	      erro12 >= errmax * .99)) {
	if (extrap)
	  ++iroff2;
	else /* if(! extrap) */
	  ++iroff1;
      }
      if (*last > 10 && all(erro12 > errmax)) // any or all?
	++iroff3;
    }
    rlist.col(maxerr-1) = area1;
    rlist.col(*last-1) = area2;
    errbnd = max(ones(ny)*epsabs, epsrel * abs(area));

    /*           test for roundoff error and eventually set error flag. */
    if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
      *ier = 2;
    if (iroff2 >= 5)
      ierro = 3;

    /* set error flag in the case that the number of subintervals equals limit. */
    if (*last == limit)
      *ier = 1;

    /*           set error flag in the case of bad integrand behaviour
		 at a point of the integration range. */

    if (R::fmax2(fabs(a1), fabs(b2)) <=
	(epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) {
      *ier = 4;
    }

    /*           append the newly-created intervals to the list. */
    if (max(error2) > max(error1)) {
      alist[maxerr] = a2;
      alist[*last] = a1;
      blist[*last] = b1;
      rlist.col(maxerr-1) = area2;
      rlist.col(*last-1) = area1;
      elist.col(maxerr-1) = error2;
      elist.col(*last-1) = error1;
      elistSum[maxerr-1] = max(error2);
      elistSum[*last-1] = max(error1);
    } else {
      alist[*last] = a2;
      blist[maxerr] = b1;
      blist[*last] = b2;
      elist.col(maxerr-1) = error1;
      elist.col(*last-1) = error2;
      elistSum[maxerr-1] = max(error1);
      elistSum[*last-1] = max(error2);
    }

    /*           call subroutine dqpsrt to maintain the descending ordering
		 in the list of error estimates and select the subinterval
		 with nrmax-th largest error estimate (to be bisected next). */

    /*L30:*/
    // Constant values: limit, last, elist
    // Updated values: ermax, maxerr, iord, nrmax
    errmaxsum = max(errmax);
    rdqpsrt(limit, last, &maxerr, &errmaxsum, elistSum, &iord[1], &nrmax);
    errmax = elist.col(maxerr-1);

    if (all(errsum <= errbnd))   goto L115;/* ***jump out of do-loop */
    if (*ier != 0)		break;
    if (*last == 2)	{ /* L80: */
      small = fabs(b - a) * .375;
      erlarg = errsum;
      ertest = errbnd;
      for (int i=0; i<ny; ++i)
	rlist2[i*52+1] = area(i);
      continue;
    }
    if (noext)		continue;

    erlarg -= erlast;
    if (fabs(b1 - a1) > small) {
      erlarg += erro12;
    }
    if (!extrap) {

      /*          test whether the interval to be bisected next is the
		  smallest interval. */

      if (fabs(blist[maxerr] - alist[maxerr]) > small) {
	continue;
      }
      extrap = TRUE;
      nrmax = 2;
    }

    if (ierro != 3 && any(erlarg > ertest)) { // any or all?

      /*           the smallest interval has the largest error.
		   before bisecting decrease the sum of the errors over the
		   larger intervals (erlarg) and perform extrapolation. */

      id = nrmax;
      jupbnd = *last;
      if (*last > limit / 2 + 2) {
	jupbnd = limit + 3 - *last;
      }
      for (k = id; k <= jupbnd; ++k) {
	maxerr = iord[nrmax];
	errmax = elist.col(maxerr-1);
	if (fabs(blist[maxerr] - alist[maxerr]) > small) {
	  goto L90;
	}
	++nrmax;
	/* L50: */
      }
    }
    /*           perform extrapolation.  L60: */

    ++numrl2;
    for (int i=0; i<ny; ++i)
      rlist2[i*52+(numrl2 - 1)] = area[i];
    for (int i=0; i<ny; ++i) {
      double resepsi = 0.0;
      double absepsi = 0.0;
      int numrl2i = numrl2;
      int nresi = nres;
      rdqelg(&numrl2i, &rlist2[i*52], &resepsi, &absepsi, &res3la[i*3], &nresi);
      reseps(i) = resepsi;
      abseps(i) = absepsi;
      if (i==ny-1) numrl2 = numrl2i;
    }
    ++ktmin;
    if (ktmin > 5 && all(abserr < errsum * .001)) {
      *ier = 5; // "roundoff error is detected in the extrapolation table"
    }
    if (any(abseps < abserr)) {
      ktmin = 0;
      abserr = abseps;
      result = reseps;
      correc = erlarg;
      ertest = max(ones(ny)*epsabs, epsrel * abs(reseps));
      if (all(abserr <= ertest)) {
	break;
      }
    }

    /*           prepare bisection of the smallest interval.  L70: */

    if (numrl2 == 1) {
      noext = TRUE;
    }
    if (*ier == 5) {
      break;
    }
    maxerr = iord[1];
    errmax = elist.col(maxerr-1);
    nrmax = 1;
    extrap = FALSE;
    small *= .5;
    erlarg = errsum;
  L90:
    ;
  }


  /* L100:	set final result and error estimate. */
  /*		------------------------------------ */

  if (all(abserr == oflow)) 	goto L115; // all or any?
  if (*ier + ierro != 0) {
    if (ierro == 3)
      abserr += correc;
    if (*ier == 0)
      *ier = 3;
    if (all(result == 0. || area == 0.)) {
      if (all(abserr > errsum)) 	goto L115;
      if (all(area == 0.)) 		goto L130;
    }
    else { /* L105:*/
      if (all(abserr / abs(result) > errsum / abs(area)))
	goto L115;
    }
  }

  /* L110: test on divergence. */
  if (!(ksgn == -1 && all(max(abs(result), abs(area)) <= defabs * .01)) &&
      all(.01 > result / area || result / area > 100. || errsum > abs(area))) {
    *ier = 5; // "the integral is probably divergent"
  }
  goto L130;

 L115:/*		compute global integral sum. */
  result = zeros(ny);
  for (k = 1; k <= *last; ++k)
    result += rlist.col(k-1);
  abserr = errsum;
 L130:
 L140:
  *neval = *last * 42 - 21;
  for (int i=0; i<ny; ++i) {
    resultp[i] = result[i]*(upper[i]-lower[i]);
    abserrp[i] = abserr[i]*(upper[i]-lower[i]);
  }
  R_Free(res3la);
  R_Free(rlist2);
  R_Free(elistSum);
  return;
} /* vrdqagse_ */

static double c_b6 = 0.;
static double c_b7 = 1.;

template<typename F>
void vrdqagie(const F f, const vec bboun, const int inf,
	      const double epsabs, const double epsrel, const int limit, const int ny,
	      double *resultp, double *abserrp, int *neval, int *ier, double *alist,
	      double *blist, double *rlistp, double *elistp, int *iord, int *last) {
  /* Local variables */
  vec area=zeros(ny), dres=zeros(ny);
  int ksgn;
  int nres;
  vec area1=zeros(ny), area2=zeros(ny), area12=zeros(ny), erro12=zeros(ny);
  int k;
  double small = 0.0;
  int ierro;
  double a1, a2, b1, b2, oflow;
  vec defab1=zeros(ny), defab2=zeros(ny);
  int ktmin, nrmax;
  double uflow;
  Rboolean noext;
  int iroff1, iroff2, iroff3;
  double *res3la = R_Calloc(3*ny, double);
  vec error1=zeros(ny), error2=zeros(ny);
  int id;
  double *rlist2 = R_Calloc(52*ny, double);
  int numrl2;
  vec correc=zeros(ny), abseps=zeros(ny), errbnd=zeros(ny), resabs=zeros(ny), erlarg=zeros(ny), defabs=zeros(ny);
  double epmach;
  int jupbnd;
  vec erlast=zeros(ny), errmax=zeros(ny), reseps=zeros(ny);
  int maxerr;
  double errmaxsum;
  Rboolean extrap;
  vec ertest=zeros(ny), errsum=zeros(ny);
  vec result = zeros(ny), abserr = zeros(ny);

  mat elist(elistp, ny, limit, false);
  mat rlist(rlistp, ny, limit, false);
  double *elistSum = R_Calloc(limit, double);

  /* ***first executable statement  dqagie */
  /* Parameter adjustments */
  --iord;
  // --elist;
  // --rlist;
  --blist;
  --alist;

  /* Function Body */
  epmach = DBL_EPSILON;

  /*           test on validity of parameters */
  /*           ----------------------------- */

  *ier = 0;
  *neval = 0;
  *last = 0;
  alist[1] = 0.;
  blist[1] = 1.;
  rlist.col(0) = zeros(ny);
  elist.col(0) = zeros(ny);
  elistSum[0] = 0.0;
  iord[1] = 0;
  if (epsabs <= 0. && (epsrel < R::fmax2(epmach * 50., 5e-29))) *ier = 6;
  if (*ier == 6) return;

  /*           first approximation to the integral */
  /*           ----------------------------------- */

  /*         determine the interval to be mapped onto (0,1).
	     if inf = 2 the integral is computed as i = i1+i2, where
	     i1 = integral of f over (-infinity,0),
	     i2 = integral of f over (0,+infinity). */

  vrdqk15i<F>(f, bboun, inf, &c_b6, &c_b7, &result, &abserr, &defabs, &resabs);

  /*           test on accuracy */

  *last = 1;
  rlist.col(0) = result;
  elist.col(0) = abserr;
  elistSum[0] = max(abserr);
  iord[1] = 1;
  dres = abs(result);
  errbnd = max(ones(ny)*epsabs, epsrel * dres);
  if (all(abserr <= epmach * 100. * defabs && abserr > errbnd))
    *ier = 2; // "roundoff error was detected"
  if (limit == 1) *ier = 1; // "maximum number of subdivisions reached"
  if (*ier != 0
      || all(abserr <= errbnd && abserr != resabs)
      || all(abserr == 0.))
    goto L130;

  /*           initialization */
  /*           -------------- */

  uflow = DBL_MIN;
  oflow = DBL_MAX;
  for (int i=0; i<ny; ++i)
    rlist2[i*52] = result[i];
  errmax = abserr;
  maxerr = 1;
  area = result;
  errsum = abserr;
  for (int i=0; i<ny; ++i)
    abserr[i] = oflow;
  nrmax = 1;
  nres = 0;
  ktmin = 0;
  numrl2 = 2;
  extrap = FALSE;
  noext = FALSE;
  ierro = 0;
  iroff1 = 0;
  iroff2 = 0;
  iroff3 = 0;
  ksgn = -1;
  if (any(dres >= (1. - epmach * 50.) * defabs)) { // any or all or do separately?
    ksgn = 1;
  }

  /*           main do-loop */
  /*           ------------ */

  for (*last = 2; *last <= limit; ++(*last)) {

    /*           bisect the subinterval with nrmax-th largest error estimate. */

    a1 = alist[maxerr];
    b1 = (alist[maxerr] + blist[maxerr]) * .5;
    a2 = b1;
    b2 = blist[maxerr];
    erlast = errmax;
    vrdqk15i<F>(f, bboun, inf, &a1, &b1, &area1, &error1, &resabs, &defab1);
    vrdqk15i<F>(f, bboun, inf, &a2, &b2, &area2, &error2, &resabs, &defab2);

    /*           improve previous approximations to integral
		 and error and test for accuracy. */

    area12 = area1 + area2;
    erro12 = error1 + error2;
    errsum = errsum + erro12 - errmax;
    area = area + area12 - rlist.col(maxerr-1);
    if (all(defab1 != error1 && defab2 != error2)) {
      if (all(abs(rlist.col(maxerr-1) - area12) <= abs(area12) * 1e-5 &&
	      erro12 >= errmax * .99)) {
	if (extrap)
	  ++iroff2;
	else /* if (! extrap) */
	  ++iroff1;
      }
      if (*last > 10 && all(erro12 > errmax))
	++iroff3;
    }

    rlist.col(maxerr-1) = area1;
    rlist.col(*last-1) = area2;
    errbnd = max(ones(ny)*epsabs, epsrel * abs(area));

    /*           test for roundoff error and eventually set error flag. */

    if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
      *ier = 2;
    if (iroff2 >= 5)
      ierro = 3;

    /*           set error flag in the case that the number of
		 subintervals equals limit. */

    if (*last == limit)
      *ier = 1;

    /*           set error flag in the case of bad integrand behaviour
		 at some points of the integration range. */

    if (R::fmax2(fabs(a1), fabs(b2)) <=
	(epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3))
      {
	*ier = 4;
      }

    /*           append the newly-created intervals to the list. */

    if (max(error2) <= max(error1)) {
      alist[*last] = a2;
      blist[maxerr] = b1;
      blist[*last] = b2;
      elist.col(maxerr-1) = error1;
      elist.col(*last-1) = error2;
      elistSum[maxerr-1] = max(error1);
      elistSum[*last-1] = max(error2);
    }
    else {
      alist[maxerr] = a2;
      alist[*last] = a1;
      blist[*last] = b1;
      rlist.col(maxerr-1) = area2;
      rlist.col(*last-1) = area1;
      elist.col(maxerr-1) = error2;
      elist.col(*last-1) = error1;
      elistSum[maxerr-1] = max(error2);
      elistSum[*last-1] = max(error1);
    }

    /*           call subroutine dqpsrt to maintain the descending ordering
		 in the list of error estimates and select the subinterval
		 with nrmax-th largest error estimate (to be bisected next). */

    errmaxsum = max(errmax);
    rdqpsrt(limit, last, &maxerr, &errmaxsum, elistSum, &iord[1], &nrmax);
    errmax = elist.col(maxerr-1);
    if (all(errsum <= errbnd)) {
      goto L115;
    }
    if (*ier != 0)	    break;
    if (*last == 2) { /* L80: */
      small = .375;
      erlarg = errsum;
      ertest = errbnd;
      for (int i=0; i<ny; ++i)
	rlist2[i*52+1] = area(i);
      continue;
    }
    if (noext) 	    continue;

    erlarg -= erlast;
    if (fabs(b1 - a1) > small) {
      erlarg += erro12;
    }
    if (!extrap) {

      /*           test whether the interval to be bisected next is the
		   smallest interval. */

      if (fabs(blist[maxerr] - alist[maxerr]) > small) {
	continue;
      }
      extrap = TRUE;
      nrmax = 2;
    }

    if (ierro != 3 && any(erlarg > ertest)) {

      /*	    the smallest interval has the largest error.
		    before bisecting decrease the sum of the errors over the
		    larger intervals (erlarg) and perform extrapolation. */

      id = nrmax;
      jupbnd = *last;
      if (*last > limit / 2 + 2) {
	jupbnd = limit + 3 - *last;
      }
      for (k = id; k <= jupbnd; ++k) {
	maxerr = iord[nrmax];
	errmax = elist.col(maxerr-1);
	if (fabs(blist[maxerr] - alist[maxerr]) > small) {
	  goto L90;
	}
	++nrmax;
	/* L50: */
      }
    }
    /*           perform extrapolation.  L60: */
    ++numrl2;
    for (int i=0; i<ny; ++i)
      rlist2[i*52+(numrl2 - 1)] = area[i];
    // rdqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
    for (int i=0; i<ny; ++i) {
      double resepsi = 0.0;
      double absepsi = 0.0;
      int numrl2i = numrl2;
      int nresi = nres;
      rdqelg(&numrl2i, &rlist2[i*52], &resepsi, &absepsi, &res3la[i*3], &nresi);
      reseps(i) = resepsi;
      abseps(i) = absepsi;
      if (i==ny-1) numrl2 = numrl2i;
    }

    ++ktmin;
    if (ktmin > 5 && all(abserr < errsum * .001)) {
      *ier = 5;
    }
    if (all(abseps >= abserr)) {
      goto L70;
    }
    ktmin = 0;
    abserr = abseps;
    result = reseps;
    correc = erlarg;
    ertest = max(ones(ny)*epsabs, epsrel * abs(reseps));
    if (all(abserr <= ertest)) {
      break;
    }

    /*            prepare bisection of the smallest interval. */

  L70:
    if (numrl2 == 1) {
      noext = TRUE;
    }
    if (*ier == 5) {
      break;
    }
    maxerr = iord[1];
    errmax = elist.col(maxerr-1);
    nrmax = 1;
    extrap = FALSE;
    small *= .5;
    erlarg = errsum;
  L90:
    ;
  }

  /* L100:     set final result and error estimate. */
  /*	     ------------------------------------ */

  if (all(abserr == oflow)) {
    goto L115;
  }
  if (*ier + ierro == 0) {
    goto L110;
  }
  if (ierro == 3) {
    abserr += correc;
  }
  if (*ier == 0) {
    *ier = 3;
  }
  if (all(result == 0. || area == 0.)) {
    if (all(abserr > errsum))
      goto L115;

    if (all(area == 0.))
      goto L130;
  }
  else { /* L105: */
    if (all(abserr / abs(result) > errsum / abs(area))) {
      goto L115;
    }
  }

  /*           test on divergence */
 L110:
  if (ksgn == -1 && all(max(abs(result), abs(area)) <= defabs * .01)) {
    goto L130;
  }
  if (all(.01 > result / area || result / area > 100. || errsum > abs(area))) {
    *ier = 6;
  }
  goto L130;

  /*           compute global integral sum. */

 L115:
  result = zeros(ny);
  for (k = 1; k <= *last; ++k)
    result += rlist.col(k-1);

  abserr = errsum;
 L130:
  *neval = *last * 30 - 15;
  if (inf == 2) {
    *neval <<= 1;
  }
  if (*ier > 2) {
    --(*ier);
  }
  for (int i=0; i<ny; ++i) {
    resultp[i] = result[i];
    abserrp[i] = abserr[i];
  }
  R_Free(res3la);
  R_Free(rlist2);
  R_Free(elistSum);
  return;
} /* vrdqagie_ */


template<typename F>
void vrdqk15i(const F f,
		     const vec bboun,
		     const int inf, double *a, double *b,
		     vec *result,
		     vec *abserr, vec *resabs, vec *resasc)
{
  /* Initialized data */

  static double wg[8] = {
    0., .129484966168869693270611432679082,
    0., .27970539148927666790146777142378,
    0., .381830050505118944950369775488975,
    0., .417959183673469387755102040816327 };
  static double xgk[8] = {
    .991455371120812639206854697526329,
    .949107912342758524526189684047851,
    .864864423359769072789712788640926,
    .741531185599394439863864773280788,
    .58608723546769113029414483825873,
    .405845151377397166906606412076961,
    .207784955007898467600689403773245, 0. };
  static double wgk[8] = {
    .02293532201052922496373200805897,
    .063092092629978553290700663189204,
    .104790010322250183839876322541518,
    .140653259715525918745189590510238,
    .16900472663926790282658342659855,
    .190350578064785409913256402421014,
    .204432940075298892414161999234649,
    .209482141084727828012999174891714 };

  /* Local variables */
  double epmach, uflow, dinf, centr, hlgth, boun;
  std::vector<vec> fv1(7), fv2(7);
  double vect[15], vec2[15];
  double tabsc1, tabsc2, absc, absc1, absc2;
  vec fval1, fval2, fc, resg, resk, fsum, reskh;
  int j, ny;

  /* ***first executable statement  dqk15i */
  epmach = DBL_EPSILON;
  uflow = DBL_MIN;
  dinf = (double) R::imin2(1, inf);
  boun = 0.0;

  centr = (*a + *b) * .5;
  hlgth = (*b - *a) * .5;
  tabsc1 = boun + dinf * (1. - centr) / centr;
  vect[0] = tabsc1;
  if (inf == 2) {
    vec2[0] = -tabsc1;
  }
  for (j = 1; j <= 7; ++j) {
    absc = hlgth * xgk[j - 1];
    absc1 = centr - absc;
    absc2 = centr + absc;
    tabsc1 = boun + dinf * (1. - absc1) / absc1;
    tabsc2 = boun + dinf * (1. - absc2) / absc2;
    vect[(j << 1) - 1] = tabsc1;
    vect[j * 2] = tabsc2;
    if (inf == 2) {
      vec2[(j << 1) - 1] = -tabsc1;
      vec2[j * 2] = -tabsc2;
    }
    /* L5: */
  }
  fval1 = eval_f<F>(f, bboun+vect[0]);
  ny = fval1.size();
  if (inf == 2) fval2 = eval_f<F>(f,bboun+vec2[0]);
  if (inf == 2) fval1 += fval2;
  fc = fval1 / centr / centr;

  /*           compute the 15-point kronrod approximation to
	       the integral, and estimate the error. */

  resg = wg[7] * fc;
  resk = wgk[7] * fc;
  *resabs = abs(resk);
  for (j = 1; j <= 7; ++j) {
    absc = hlgth * xgk[j - 1];
    absc1 = centr - absc;
    absc2 = centr + absc;
    tabsc1 = boun + dinf * (1. - absc1) / absc1;
    tabsc2 = boun + dinf * (1. - absc2) / absc2;
    fval1 = eval_f<F>(f,bboun+vect[(j << 1) - 1]);
    fval2 = eval_f<F>(f,bboun+vect[j * 2]);
    if (inf == 2) {
      fval1 += eval_f<F>(f,bboun+vec2[(j << 1) - 1]);
      fval2 += eval_f<F>(f,bboun+vec2[j * 2]);
    }
    fval1 = fval1 / absc1 / absc1;
    fval2 = fval2 / absc2 / absc2;
    fv1[j - 1] = fval1;
    fv2[j - 1] = fval2;
    fsum = fval1 + fval2;
    resg += wg[j - 1] * fsum;
    resk += wgk[j - 1] * fsum;
    *resabs += wgk[j - 1] * (abs(fval1) + abs(fval2));
    /* L10: */
  }
  reskh = resk * .5;
  *resasc = wgk[7] * abs(fc - reskh);
  for (j = 1; j <= 7; ++j) {
    *resasc += wgk[j - 1] * (abs(fv1[j - 1] - reskh) +
			     abs(fv2[j - 1] - reskh));
    /* L20: */
  }
  *result = resk * hlgth;
  *resasc *= hlgth;
  *resabs *= hlgth;
  *abserr = abs((resk - resg) * hlgth);
  for (int i=0; i<ny; ++i) {
    if ((*resasc)[i] != 0. && (*abserr)[i] != 0.) {
      (*abserr)[i] = (*resasc)[i] * R::fmin2(1.0, pow((*abserr)[i] * 200. / (*resasc)[i], 1.5));
    }
    if ((*resabs)[i] > uflow / (epmach * 50.)) {
      (*abserr)[i] = R::fmax2(epmach * 50. * (*resabs)[i], (*abserr)[i]);
    }
  }
  return;
} /* vrdqk15i_ */

template<typename F>
void vrdqk21(const F f, const vec lower, const vec upper, const double a,
	     const double b, vec *result,
	     vec *abserr, vec *resabs, vec *resasc)
{
  /* Initialized data */

  static double wg[5] = { .066671344308688137593568809893332,
    .149451349150580593145776339657697,
    .219086362515982043995534934228163,
    .269266719309996355091226921569469,
    .295524224714752870173892994651338 };
  static double xgk[11] = { .995657163025808080735527280689003,
    .973906528517171720077964012084452,
    .930157491355708226001207180059508,
    .865063366688984510732096688423493,
    .780817726586416897063717578345042,
    .679409568299024406234327365114874,
    .562757134668604683339000099272694,
    .433395394129247190799265943165784,
    .294392862701460198131126603103866,
    .14887433898163121088482600112972,0. };
  static double wgk[11] = { .011694638867371874278064396062192,
    .03255816230796472747881897245939,
    .05475589657435199603138130024458,
    .07503967481091995276704314091619,
    .093125454583697605535065465083366,
    .109387158802297641899210590325805,
    .123491976262065851077958109831074,
    .134709217311473325928054001771707,
    .142775938577060080797094273138717,
    .147739104901338491374841515972068,
    .149445554002916905664936468389821 };

  /* Local variables */
  std::vector<vec> fv1(10), fv2(10);
  double absc, vect[21];
  vec resg, resk, fsum, fval1, fval2;
  vec fc, reskh, delta;
  double hlgth, centr, uflow;
  double epmach, dhlgth;
  int j, jtw, jtwm1, ny;

  epmach = DBL_EPSILON;
  uflow = DBL_MIN;

  delta = upper-lower;
  centr = (a + b) * .5;
  hlgth = (b - a) * .5;
  dhlgth = fabs(hlgth);

  /*           compute the 21-point kronrod approximation to
	       the integral, and estimate the absolute error. */

  vect[0] = centr;
  for (j = 1; j <= 5; ++j) {
    jtw = j << 1;
    absc = hlgth * xgk[jtw - 1];
    vect[(j << 1) - 1] = centr - absc;
    /* L5: */
    vect[j * 2] = centr + absc;
  }
  for (j = 1; j <= 5; ++j) {
    jtwm1 = (j << 1) - 1;
    absc = hlgth * xgk[jtwm1 - 1];
    vect[(j << 1) + 9] = centr - absc;
    vect[(j << 1) + 10] = centr + absc;
  }
  fc = eval_f<F>(f,lower+delta*vect[0]);
  ny = fc.size();
  resg = zeros(ny);
  resk = wgk[10] * fc;
  *resabs = abs(resk);
  for (j = 1; j <= 5; ++j) {
    jtw = j << 1;
    absc = hlgth * xgk[jtw - 1];
    fval1 = eval_f<F>(f,lower+delta*vect[(j << 1) - 1]);
    fval2 = eval_f<F>(f,lower+delta*vect[j * 2]);
    fv1[jtw - 1] = fval1;
    fv2[jtw - 1] = fval2;
    fsum = fval1 + fval2;
    resg += wg[j - 1] * fsum;
    resk += wgk[jtw - 1] * fsum;
    *resabs += wgk[jtw - 1] * (abs(fval1) + abs(fval2));
    /* L10: */
  }
  for (j = 1; j <= 5; ++j) {
    jtwm1 = (j << 1) - 1;
    absc = hlgth * xgk[jtwm1 - 1];
    fval1 = eval_f<F>(f,lower+delta*vect[(j << 1) + 9]);
    fval2 = eval_f<F>(f,lower+delta*vect[(j << 1) + 10]);
    fv1[jtwm1 - 1] = fval1;
    fv2[jtwm1 - 1] = fval2;
    fsum = fval1 + fval2;
    resk += wgk[jtwm1 - 1] * fsum;
    *resabs += wgk[jtwm1 - 1] * (abs(fval1) + abs(fval2));
    /* L15: */
  }
  reskh = resk * .5;
  *resasc = wgk[10] * abs(fc - reskh);
  for (j = 1; j <= 10; ++j) {
    *resasc += wgk[j - 1] * (abs(fv1[j - 1] - reskh) +
			     abs(fv2[j - 1] - reskh));
    /* L20: */
  }
  *result = resk * hlgth;
  *resabs *= dhlgth;
  *resasc *= dhlgth;
  *abserr = abs((resk - resg) * hlgth);
  for (int i=0; i<ny; ++i) {
    if ((*resasc)[i] != 0. && (*abserr)[i] != 0.) {
      (*abserr)[i] = (*resasc)[i] * R::fmin2(1.0, pow((*abserr)[i] * 200. / (*resasc)[i], 1.5));
    }
    if ((*resabs)[i] > uflow / (epmach * 50.)) {
      (*abserr)[i] = R::fmax2(epmach * 50. * (*resabs)[i], (*abserr)[i]);
    }
  }
  return;
} /* vrdqk21_ */

// Constant values: limit, last, elist
// Updated values: ermax, maxerr, iord, nrmax
static void rdqpsrt(const int limit, int *last, int *maxerr,
		    double *ermax, double *elist, int *iord, int *nrmax)
{
  /* Local variables */
  int i, j, k, ido, jbnd, isucc, jupbn;
  double errmin, errmax;

  /* Parameter adjustments */
  --iord;
  --elist;

  /* Function Body */

  /*           check whether the list contains more than
	       two error estimates. */
  if (*last <= 2) {
    iord[1] = 1;
    iord[2] = 2;
    goto Last;
  }
  /*           this part of the routine is only executed if, due to a
	       difficult integrand, subdivision increased the error
	       estimate. in the normal case the insert procedure should
	       start after the nrmax-th largest error estimate. */

  errmax = elist[*maxerr];
  if (*nrmax > 1) {
    ido = *nrmax - 1;
    for (i = 1; i <= ido; ++i) {
      isucc = iord[*nrmax - 1];
      if (errmax <= elist[isucc])
	break; /* out of for-loop */
      iord[*nrmax] = isucc;
      --(*nrmax);
      /* L20: */
    }
  }

  /*L30:       compute the number of elements in the list to be maintained
    in descending order. this number depends on the number of
    subdivisions still allowed. */
  if (*last > limit / 2 + 2)
    jupbn = limit + 3 - *last;
  else
    jupbn = *last;

  errmin = elist[*last];

  /*           insert errmax by traversing the list top-down,
	       starting comparison from the element elist(iord(nrmax+1)). */

  jbnd = jupbn - 1;
  for (i = *nrmax + 1; i <= jbnd; ++i) {
    isucc = iord[i];
    if (errmax >= elist[isucc]) {/* ***jump out of do-loop */
      /* L60: insert errmin by traversing the list bottom-up. */
      iord[i - 1] = *maxerr;
      for (j = i, k = jbnd; j <= jbnd; j++, k--) {
	isucc = iord[k];
	if (errmin < elist[isucc]) {
	  /* goto L80; ***jump out of do-loop */
	  iord[k + 1] = *last;
	  goto Last;
	}
	iord[k + 1] = isucc;
      }
      iord[i] = *last;
      goto Last;
    }
    iord[i - 1] = isucc;
  }

  iord[jbnd] = *maxerr;
  iord[jupbn] = *last;

 Last:/* set maxerr and ermax. */

  *maxerr = iord[*nrmax];
  *ermax = elist[*maxerr];
  return;
} /* rdqpsrt_ */

// Constant values:
// Updated values: n, epstab, result, abserr, res3la, nres
static void rdqelg(int *n, double *epstab, double *
		   result, double *abserr, double *res3la, int *nres)
{
  /* Local variables */
  int i__, indx, ib, ib2, ie, k1, k2, k3, num, newelm, limexp;
  double delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epmach, epsinf;
  double oflow, ss, res;
  double errA, err1, err2, err3, tol1, tol2, tol3;

  /* Parameter adjustments */
  --res3la;
  --epstab;

  /* Function Body */
  epmach = DBL_EPSILON;
  oflow = DBL_MAX;
  ++(*nres);
  *abserr = oflow;
  *result = epstab[*n];
  if (*n < 3) {
    goto L100;
  }
  limexp = 50;
  epstab[*n + 2] = epstab[*n];
  newelm = (*n - 1) / 2;
  epstab[*n] = oflow;
  num = *n;
  k1 = *n;
  for (i__ = 1; i__ <= newelm; ++i__) {
    k2 = k1 - 1;
    k3 = k1 - 2;
    res = epstab[k1 + 2];
    e0 = epstab[k3];
    e1 = epstab[k2];
    e2 = res;
    e1abs = fabs(e1);
    delta2 = e2 - e1;
    err2 = fabs(delta2);
    tol2 = R::fmax2(fabs(e2), e1abs) * epmach;
    delta3 = e1 - e0;
    err3 = fabs(delta3);
    tol3 = R::fmax2(e1abs, fabs(e0)) * epmach;
    if (err2 <= tol2 && err3 <= tol3) {
      /*           if e0, e1 and e2 are equal to within machine
		   accuracy, convergence is assumed. */
      *result = res;/*		result = e2 */
      *abserr = err2 + err3;/*	abserr = fabs(e1-e0)+fabs(e2-e1) */

      goto L100;	/* ***jump out of do-loop */
    }

    e3 = epstab[k1];
    epstab[k1] = e1;
    delta1 = e1 - e3;
    err1 = fabs(delta1);
    tol1 = R::fmax2(e1abs, fabs(e3)) * epmach;

    /*           if two elements are very close to each other, omit
		 a part of the table by adjusting the value of n */

    if (err1 > tol1 && err2 > tol2 && err3 > tol3) {
      ss = 1. / delta1 + 1. / delta2 - 1. / delta3;
      epsinf = fabs(ss * e1);

      /*           test to detect irregular behaviour in the table, and
		   eventually omit a part of the table adjusting the value of n. */

      if (epsinf > 1e-4) {
	goto L30;
      }
    }

    *n = i__ + i__ - 1;
    goto L50;/* ***jump out of do-loop */


  L30:/* compute a new element and eventually adjust the value of result. */

    res = e1 + 1. / ss;
    epstab[k1] = res;
    k1 += -2;
    errA = err2 + fabs(res - e2) + err3;
    if (errA <= *abserr) {
      *abserr = errA;
      *result = res;
    }
  }

  /*           shift the table. */

 L50:
  if (*n == limexp) {
    *n = (limexp / 2 << 1) - 1;
  }

  if (num / 2 << 1 == num) ib = 2; else ib = 1;
  ie = newelm + 1;
  for (i__ = 1; i__ <= ie; ++i__) {
    ib2 = ib + 2;
    epstab[ib] = epstab[ib2];
    ib = ib2;
  }
  if (num != *n) {
    indx = num - *n + 1;
    for (i__ = 1; i__ <= *n; ++i__) {
      epstab[i__] = epstab[indx];
      ++indx;
    }
  }
  /*L80:*/
  if (*nres >= 4) {
    /* L90: */
    *abserr = fabs(*result - res3la[3]) +
      fabs(*result - res3la[2]) +
      fabs(*result - res3la[1]);
    res3la[1] = res3la[2];
    res3la[2] = res3la[3];
    res3la[3] = *result;
  } else {
    res3la[*nres] = *result;
    *abserr = oflow;
  }

 L100:/* compute error estimate */
  *abserr = R::fmax2(*abserr, epmach * 5. * fabs(*result));
  return;
} /* rdqelg_ */


template<typename F>
Rcpp::List vdqags(const F f, const vec a, const vec b,
		  const double epsrel, const double epsabs, const int limit,
		  const int ny) {
  using namespace Rcpp;
  double *result = R_Calloc(ny, double);
  double *abserr = R_Calloc(ny, double);
  int lenw, ier, neval, last;
  lenw = 2*limit*ny + 2*limit;
  int *iwork = R_Calloc(limit, int);
  double *work = R_Calloc(lenw, double);
  vRdqags<F>(f, a, b, epsabs, epsrel, ny,
	     result, abserr, &neval, &ier, limit, &lenw, &last, iwork, work);
  vec result2(result,ny);
  vec abserr2(abserr,ny);
  R_Free(result);
  R_Free(abserr);
  R_Free(iwork);
  R_Free(work);
  return List::create(_("value") = result2,
		      _("abs.err") = abserr2,
		      _("subdivisions") = last,
		      _("ierr") = ier);
}

template<typename F>
Rcpp::List vdqagi(const F f, const vec bound, const int inf,
		  const double epsrel, const double epsabs, const int limit, const int ny) {
  using namespace Rcpp;
  double *result = R_Calloc(ny, double);
  double *abserr = R_Calloc(ny, double);
  int lenw, ier, neval, last;
  lenw = 2*limit*ny + 2*limit;
  int *iwork = R_Calloc(limit, int);
  double *work = R_Calloc(lenw, double);
  vRdqagi(f, bound, inf, epsabs, epsrel, limit, ny,
	  result, abserr, &neval, &ier, &lenw, &last, iwork, work);
  vec resultnv(result, ny);
  vec abserrnv(abserr, ny);
  R_Free(result);
  R_Free(abserr);
  R_Free(iwork);
  R_Free(work);
  return List::create(_("value") = resultnv,
		      _("abs.err") = abserrnv,
		      _("subdivisions") = last,
		      _("ierr") = ier);
}

// [[Rcpp::export]]
Rcpp::List vdqagsRcpp(const Rcpp::Function f, const arma::vec a, const arma::vec b,
		      const double epsrel, const double epsabs, const int limit,
		      const int ny) {
// RcppExport SEXP vdqagsRcpp(SEXP _f, SEXP _a, SEXP _b,
// 			   SEXP _epsrel, SEXP _epsabs, SEXP _limit,
// 			   SEXP _ny) {
  using namespace Rcpp;
  // Function f = as<Function>(_f);
  // double a = as<double>(_a);
  // double b = as<double>(_b);
  // double epsrel = as<double>(_epsrel);
  // double epsabs = as<double>(_epsabs);
  // int limit = as<int>(_limit);
  // int ny = as<int>(_ny);
  return vdqags(f, a, b, epsrel, epsabs, limit, ny);
}

// [[Rcpp::export]]
Rcpp::List vdqagiRcpp(const Rcpp::Function f, const arma::vec bound, const int inf,
		      const double epsrel, const double epsabs, const int limit, const int ny) {
// RcppExport SEXP vdqagiRcpp(SEXP _f, SEXP _bound, SEXP _inf,
// 			   SEXP _epsrel, SEXP _epsabs, SEXP _limit, SEXP _ny) {
  using namespace Rcpp;
  // Function f = as<Function>(_f);
  // double bound = as<double>(_bound);
  // int inf = as<int>(_inf);
  // double epsrel = as<double>(_epsrel);
  // double epsabs = as<double>(_epsabs);
  // int limit = as<int>(_limit);
  // int ny = as<int>(_ny);
  return vdqagi(f, bound, inf, epsrel, epsabs, limit, ny);
}

// [[Rcpp::export]]
Rcpp::List vrdqk21Rcpp(const Rcpp::Function f, const arma::vec lower, const arma::vec upper,
		       const double a, const double b) {
// RcppExport SEXP vrdqk21Rcpp(SEXP _f, SEXP _a, SEXP _b) {
  using namespace Rcpp;
  // Function f = as<Function>(_f);
  // double a = as<double>(_a);
  // double b = as<double>(_b);
  vec result, abserr, resabs, resasc;
  vrdqk21(f, lower, upper, a, b, &result, &abserr, &resabs, &resasc);
  return List::create(_("value") = result,
		      _("abs.err") = abserr);
}

// [[Rcpp::export]]
Rcpp::List vrdqk15Rcpp(const Rcpp::Function f, const arma::vec boun, const int inf,
		       double a, double b) {
// RcppExport SEXP vrdqk15Rcpp(SEXP _f, SEXP _boun, SEXP _inf, SEXP _a, SEXP _b) {
  using namespace Rcpp;
  // Function f = as<Function>(_f);
  // double boun = as<double>(_boun);
  // int inf = as<int>(_inf);
  // double a = as<double>(_a);
  // double b = as<double>(_b);
  vec result, abserr, resabs, resasc;
  vrdqk15i(f, boun, inf, &a, &b, &result, &abserr, &resabs, &resasc);
  return List::create(_("value") = result,
		      _("abserr") = abserr);
}

vec test_f(const arma::vec a) {
  vec out{exp(a[0]), exp(2*a[1]), log(a[2])};
  return out;
}

// [[Rcpp::export]]
Rcpp::List test_vdqags() {
  vec lower{0.0,0.0,0.0}, upper{1.0, 1.0, 1.0};
  return vdqags(test_f, lower, upper, 1.0e-8, 1.0e-8, 50, 3);
}

vec test_f2(const arma::vec a) {
  vec out{exp(a[0]), exp(2*a[1])};
  return out;
}

// [[Rcpp::export]]
Rcpp::List test_vdqagi() {
  vec bound{0.0,0.0};
  return vdqagi(test_f2, bound, -1, 1.0e-8, 1.0e-8, 50, 2);
}
