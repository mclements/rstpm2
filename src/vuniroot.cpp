#include <cfloat>
#include <Rcpp.h>

//[[Rcpp::export]]
Rcpp::List vunirootRcpp(Rcpp::Function f,
			Rcpp::NumericVector lower,
			Rcpp::NumericVector upper,
			Rcpp::NumericVector fa,
			Rcpp::NumericVector fb,
			int numiter,
			double tol) {
  using namespace Rcpp;
    int size = lower.size();
    NumericVector a(clone(lower)), b(clone(upper)), c(clone(a)), Tol(size,0.0);
    NumericVector fc(clone(fa));
    LogicalVector converged(size,false);
    IntegerVector ns(size,-1);
    int i;
    /* First test if we have found a root at an endpoint */
    for(i=0; i<size; i++) {
	if (fa[i]==0.0) {
	  converged[i]=true;
	  b[i]=a[i];
	}
	if (fb[i]==0.0) {
	  converged[i]=true;
	}
    }
    for (int n = 1; n<=numiter; n++) {
      R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
      for(i=0; i<size; i++) {
	if (!converged[i]) {
	  double prev_step = b[i]-a[i];		/* Distance from the last but one
						   to the last approximation	*/
	  double tol_act;			/* Actual tolerance		*/
	  double p;			/* Interpolation step is calcu- */
	  double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	  double new_step;		/* Step at this iteration	*/
	  if( fabs(fc[i]) < fabs(fb[i]) )
	    {				/* Swap data for b to be the	*/
	      a[i] = b[i];  b[i] = c[i];  c[i] = a[i];	/* best approximation		*/
	      fa[i]=fb[i];  fb[i]=fc[i];  fc[i]=fa[i];
	    }
	  tol_act = 2*DBL_EPSILON*fabs(b[i]) + tol/2;
	  new_step = (c[i]-b[i])/2.0;
	  if( fabs(new_step) <= tol_act || fb[i] == (double)0 )
	    {
	      // *Maxit -= maxit;
	      Tol[i] = fabs(c[i]-b[i]);
	      converged[i] = true;			/* Acceptable approx. is found	*/
	      ns[i] = n;
	    } else {
	    /* Decide if the interpolation can be tried	*/
	    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
		&& fabs(fa[i]) > fabs(fb[i]) ) {	/* and was in true direction,
							 * Interpolation may be tried	*/
	      double t1,cb,t2;
	      cb = c[i]-b[i];
	      if( a[i]==c[i] ) {		/* If we have only two distinct	*/
		/* points linear interpolation	*/
		t1 = fb[i]/fa[i];		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	      }
	      else {			/* Quadric inverse interpolation*/
		q = fa[i]/fc[i];  t1 = fb[i]/fc[i];	 t2 = fb[i]/fa[i];
		p = t2 * ( cb*q*(q-t1) - (b[i]-a[i])*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	      }
	      if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	      else			/* and assign possible minus to	*/
		p = -p;			/* q				*/
	      if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		  && p < fabs(prev_step*q/2) )	/* and isnt too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	    }
	    if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	      if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	      else
		new_step = -tol_act;
	    }
	    a[i] = b[i];	fa[i] = fb[i];			/* Save the previous approx. */
	    b[i] += new_step;	// fb = (*f)(b, info);	/* Do step to a new approxim. */
	  }
	}
      }
      if (is_true(all(converged))) break;
      fb = f(b);
      for(i=0; i<size; i++) {
	if( (fb[i] > 0 && fc[i] > 0) || (fb[i] < 0 && fc[i] < 0) ) {
	  /* Adjust c for it to have a sign opposite to that of b */
	  c[i] = a[i];  fc[i] = fa[i];
	}
      }
    }
    if (is_false(all(converged)))
      for (i=0; i<size; i++) 
	if (!converged[i]) {
	  Tol[i]=fabs(c[i]-b[i]);
	  ns[i] = -1;
	}
    return Rcpp::List::create(_("root")  = b, _("iter") = ns, _("tol")=Tol);
  }


//[[Rcpp::export]]
Rcpp::NumericVector voptimizeRcpp(Rcpp::Function f,
				    Rcpp::NumericVector ax,
				    Rcpp::NumericVector bx,
				    double tol) {
    using Rcpp::NumericVector;
    using Rcpp::clone;
    using Rcpp::as;
    const int size = ax.size();
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    NumericVector a = clone(ax), b=clone(bx);
    NumericVector w, x, fx, fu, fv, fw;
    double p, q, r;
    double eps, tol3;
    NumericVector d(size, 0.0), e(size, 0.0), u(size, 0.0), xm(size,0.0),
      t2(size,0.0), tol1(size,0.0), v(size,0.0);

    /*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    for (int i=0; i<size; ++i)
      v[i] = a[i] + c * (b[i] - a[i]);
    w = clone(v);
    x = clone(v);

    fx = f(x);
    fv = clone(fx);
    fw = clone(fx);
    tol3 = tol / 3.;

    /*  main loop starts here ----------------------------------- */

    for(;;) {
      
      R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
      
      for (int i=0; i<size; ++i) {
	xm[i] = (a[i] + b[i]) * .5;
	tol1[i] = eps * fabs(x[i]) + tol3;
	t2[i] = tol1[i] * 2.;
      }

      /* check stopping criterion -- do any fail? */
      bool any_failures = false;
      for (int i=0; i<size; ++i) {
	if (fabs(x[i] - xm[i]) > t2[i] - (b[i] - a[i]) * .5) {
	  any_failures = true;
	  break;
	}
      }
      if (!any_failures) break;
	
      for (int i=0; i<size; ++i) {
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e[i]) > tol1[i]) { /* fit parabola */

	  r = (x[i] - w[i]) * (fx[i] - fv[i]);
	  q = (x[i] - v[i]) * (fx[i] - fw[i]);
	  p = (x[i] - v[i]) * q - (x[i] - w[i]) * r;
	  q = (q - r) * 2.;
	  if (q > 0.) p = -p; else q = -q;
	  r = e[i];
	  e[i] = d[i];
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a[i] - x[i]) || p >= q * (b[i] - x[i])) { /* a golden-section step */

	  if (x[i] < xm[i]) e[i] = b[i] - x[i]; else e[i] = a[i] - x[i];
	  d[i] = c * e[i];
	}
	else { /* a parabolic-interpolation step */

	  d[i] = p / q;
	  u[i] = x[i] + d[i];

	  /* f must not be evaluated too close to ax or bx */

	  if (u[i] - a[i] < t2[i] || b[i] - u[i] < t2[i]) {
	    d[i] = tol1[i];
	    if (x[i] >= xm[i]) d[i] = -d[i];
	  }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d[i]) >= tol1[i])
	  u[i] = x[i] + d[i];
	else if (d[i] > 0.)
	  u[i] = x[i] + tol1[i];
	else
	  u[i] = x[i] - tol1[i];
      }

      fu = f(u);

      /*  update  a, b, v, w, and x */

      for (int i=0; i<size; ++i) {
	if (fu[i] <= fx[i]) {
	  if (u[i] < x[i]) b[i] = x[i]; else a[i] = x[i];
	  v[i] = w[i];    w[i] = x[i];   x[i] = u[i];
	  fv[i] = fw[i]; fw[i] = fx[i]; fx[i] = fu[i];
	} else {
	  if (u[i] < x[i]) a[i] = u[i]; else b[i] = u[i];
	  if (fu[i] <= fw[i] || w[i] == x[i]) {
	    v[i] = w[i]; fv[i] = fw[i];
	    w[i] = u[i]; fw[i] = fu[i];
	  } else if (fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i]) {
	    v[i] = u[i]; fv[i] = fu[i];
	  }
	}
      }
    }
    /* end of main loop */
    
    return x;
  }
  
