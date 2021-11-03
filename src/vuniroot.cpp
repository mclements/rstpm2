#include <cfloat>
#include <Rcpp.h>

namespace rstpm2 {

  using namespace R;
  using namespace Rcpp;
  
  // List vunirootRcpp(Function f, NumericVector lower, NumericVector upper, int numiter, double tol) {
  RcppExport SEXP vunirootRcpp(SEXP __f, SEXP __lower, SEXP __upper, SEXP __fa, SEXP __fb, SEXP __numiter, SEXP __tol) {
    Rcpp::Function f = as<Rcpp::Function>(__f);
    NumericVector lower = clone(as<NumericVector>(__lower));
    NumericVector upper = clone(as<NumericVector>(__upper));
    NumericVector fa = clone(as<NumericVector>(__fa));
    NumericVector fb = clone(as<NumericVector>(__fb));
    int numiter = as<int>(__numiter);
    double tol = as<double>(__tol);
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
    return wrap(List::create(_("root")  = b, _("iter") = ns, _("tol")=Tol));
  }
  
} // namespace
