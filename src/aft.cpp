#include <RcppArmadillo.h>
#include "c_optim.h"

namespace rstpm2 {
  using namespace arma;
  using namespace Rcpp;

  mat qr_q(const mat& X, double tol = 1E-12) {
    Rcpp::NumericMatrix nmX = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(X));
    Rcpp::NumericMatrix nmQ = qr_q(nmX, tol);
    return Rcpp::as<mat>(Rcpp::wrap(nmQ));
  }
      
  class SplineBasis {
  public:
    int order,			/* order of the spline */
      ordm1,			/* order - 1 (3 for cubic splines) */
      nknots,			/* number of knots */
      curs,			/* current position in knots vector */
      boundary,		/* must have knots[(curs) <= x < knots(curs+1) */
      ncoef;			/* number of coefficients */
    /* except for the boundary case */
    vec ldel;  	/* differences from knots on the left */
    vec rdel;	/* differences from knots on the right */
    vec knots;	/* knot vector */
    vec coeff;	/* coefficients */
    vec a;		/* scratch array */
    SplineBasis(int order = 4) : order(order) {
      ordm1 = order - 1;
      rdel = vec(ordm1);
      ldel = vec(ordm1);
      a = vec(order);
    }
    SplineBasis(vec knots, int order = 4) : order(order), knots(knots) {
      ordm1 = order - 1;
      nknots = knots.size();
      ncoef = nknots - order;
      rdel = vec(ordm1);
      ldel = vec(ordm1);
      a = vec(order);
    }
    int set_cursor(double x)
    {
      int i;
      /* don't assume x's are sorted */
      curs = -1; /* Wall */
      boundary = 0;
      for (i = 0; i < nknots; i++) {
	if (knots(i) >= x) curs = i;
	if (knots(i) > x) break;
      }
      if (curs > nknots - order) {
	int lastLegit = nknots - order;
	if (x == knots(lastLegit)) {
	  boundary = 1; curs = lastLegit;
	}
      }
      return curs;
    }
    void
    diff_table(double x, int ndiff)
    {
      int i;
      for (i = 0; i < ndiff; i++) {
	rdel(i) = knots(curs + i) - x;
	ldel(i) = x - knots(curs - (i + 1));
      }
    }
    double slow_evaluate(double x, int nder)
    {
      int ti = curs, 
	lpt, apt, rpt, inner,
	outer = ordm1;
      if (boundary && nder == ordm1) { /* value is arbitrary */
	return double(0);
      }
      while(nder--) {  // FIXME: divides by zero
	for(inner = outer, apt = 0, lpt = ti - outer; inner--; apt++, lpt++)
	  a(apt) = double(outer) * (a(apt + 1) - a(apt))/(knots(lpt + outer) - knots(lpt));
	outer--;
      }
      diff_table(x, outer);
      while(outer--)
	for(apt = 0, lpt = outer, rpt = 0, inner = outer + 1;
	    inner--; lpt--, rpt++, apt++)
	  // FIXME: divides by zero
	  a(apt) = (a(apt + 1) * ldel(lpt) + a(apt) * rdel(rpt))/(rdel(rpt) + ldel(lpt));
      return a(0);
    }
    /* fast evaluation of basis functions */
    vec basis_funcs(double x)
    {
      vec b(order);
      diff_table(x, ordm1);
      b(0) = double(1);
      for (size_t j = 1; j <= (size_t)ordm1; j++) {
	double saved = double(0);
	for (size_t r = 0; r < j; r++) { // do not divide by zero
	  double den = rdel(r) + ldel(j - 1 - r);
	  if(den != double(0)) {
	    double term = b(r)/den;
	    b(r) = saved + rdel(r) * term;
	    saved = ldel(j - 1 - r) * term;
	  } else {
	    if(r != double(0) || rdel(r) != double(0))
	      b(r) = saved;
	    saved = double(0);
	  }
	}
	b(j) = saved;
      }
      return b;
    }
    vec eval(double x, int ders=0) {
      vec val(ncoef);
      val = zeros(ncoef);
      set_cursor(x);
      int io = curs - order;
      if (io < 0 || io > nknots) {
	for (size_t j = 0; j < (size_t)order; j++) {
	  val(j+io) = double(0); // R_NaN;
	}
      } else if (ders > 0) { /* slow method for derivatives */
	for(size_t i = 0; i < (size_t)order; i++) {
	  for(size_t j = 0; j < (size_t)order; j++) a(j) = double(0);
	  a(i) = double(1);
	  val(i+io) = slow_evaluate(x, ders);
	}
      } else { 		/* fast method for value */
	vec valtmp = basis_funcs(x);
	for (size_t i=0; i<valtmp.size(); i++)
	  val(i+io)=valtmp(i);
      }
      return val;
    }
    mat basis(vec x, int ders=0) {
      mat mat(x.size(), ncoef);
      for (size_t i=0; i<x.size(); i++) {
	vec vec = eval(x(i), ders);
	for  (size_t j=0; j<vec.size(); j++)
	  mat(i,j)=vec(j);
      }
      return mat;
    }
  };
  
  class bs : public SplineBasis {
  public:
    vec boundary_knots, interior_knots;
    int intercept, df;
    bs() {} // default constructor
    bs(vec boundary_knots, vec interior_knots, int intercept = 0) :
      SplineBasis(4), boundary_knots(boundary_knots), interior_knots(interior_knots),
      intercept(intercept) {
      df = intercept + 3 + interior_knots.size();
      this->nknots = interior_knots.size()+8;
      this->ncoef = this->nknots - this->order;
      this->knots = vec(this->nknots);
      for(size_t i=0; i<4;i++) {
	this->knots(i)=boundary_knots(0);
	this->knots(this->nknots-i-1)=boundary_knots(1);
      }
      if (interior_knots.size() > 0) 
	for(size_t i=0; i<interior_knots.size();i++) 
	  this->knots(i+4)=interior_knots(i);
    }      
    vec eval(double x, int ders=0) {
      vec vec;
      if (x<boundary_knots(0)) {
	double k_pivot = double(0.75)*boundary_knots(0)+double(0.25)*interior_knots(0);
	double delta = x - k_pivot;
	vec = bs::eval(k_pivot,0) +
	  bs::eval(k_pivot,1)*delta +
	  bs::eval(k_pivot,2)*delta*delta/2. +
	  bs::eval(k_pivot,3)*delta*delta*delta/6.;
      }
      else if (x>boundary_knots(1)) {
	double k_pivot = double(0.75)*boundary_knots(1)+double(0.25)*interior_knots(interior_knots.size()-1);
	double delta = x - k_pivot;
	vec = bs::eval(k_pivot,0) +
	  bs::eval(k_pivot,1)*delta +
	  bs::eval(k_pivot,2)*delta*delta/2. +
	  bs::eval(k_pivot,3)*delta*delta*delta/6.;
      }
      else  {
	vec = SplineBasis::eval(x, ders).subvec(1-intercept,df-1);
      }
      return vec;
    }
    mat basis(vec x, int ders=0) {
      mat mat(x.size(), df);
      for (size_t i=0; i<x.size(); i++) {
	vec vec = bs::eval(x(i), ders);
	for  (size_t j=0; j<vec.size(); j++)
	  mat(i,j)=vec(j);
      }
      return mat;
    }
  };

  class ns : public bs {
  public:
    vec tl0, tl1, tr0, tr1;
    mat q_matrix;
    ns() {} // default constructor
    // ns(vec boundary_knots, vec interior_knots, int intercept=0) :
    //   bs(boundary_knots, interior_knots, intercept) {
    //   // calculate the Q matrix
    //   mat const_basis = bs::basis(boundary_knots, 2);
    //   mat qd = qr_q(const_basis.t());
    //   mat qsub(qd.n_rows, qd.n_cols-2);
    //   for (size_t i=0; i<qsub.n_rows; i++)
    // 	for (size_t j=0; j<qsub.n_cols; j++)
    // 	  qsub(i,j) = qd(i,j+2);
    //   q_matrix = qsub.t();
    //   tl0 = q_matrix * bs::eval(boundary_knots(0), 0);
    //   tl1 = q_matrix * bs::eval(boundary_knots(0), 1);
    //   tr0 = q_matrix * bs::eval(boundary_knots(1), 0);
    //   tr1 = q_matrix * bs::eval(boundary_knots(1), 1);
    // }
    ns(vec boundary_knots, vec interior_knots, mat _q_matrix, int intercept=0) :
      bs(boundary_knots, interior_knots, intercept), q_matrix(_q_matrix) {
      tl0 = q_matrix * bs::eval(boundary_knots(0), 0);
      tl1 = q_matrix * bs::eval(boundary_knots(0), 1);
      tr0 = q_matrix * bs::eval(boundary_knots(1), 0);
      tr1 = q_matrix * bs::eval(boundary_knots(1), 1);
    }
    vec eval(double x, int der) {
      if(x < this->boundary_knots(0)) {
	if (der==0) 
	  return tl0 + (x - this->boundary_knots(0))*tl1;
	else if (der==1)
	  return tl1;
	else return tl1*double(0);
      } else if (x > this->boundary_knots(1)) {
	if (der==0)
	  return tr0 + (x - this->boundary_knots(1))*tr1;
	else if (der==1)
	  return tr1;
	else return tr1*double(0);
      }
      else return q_matrix * bs::eval(x,der);
    }
    mat basis(vec x, int ders=0) {
      mat mat(x.size(), this->df-2);
      for (size_t i=0; i<x.size(); i++) {
	vec vec = ns::eval(x(i), ders);
	for  (size_t j=0; j<vec.size(); j++)
	  mat(i,j)=vec(j);
      }
      return mat;
    }
  }; // class ns

  class aft {
  public:
    List args;
    vec init;
    mat X, X0;
    mat XD;
    vec event;
    vec time, time0;
    vec boundaryKnots;
    vec interiorKnots;
    mat q_matrix;
    ns s;
    bool delayed;
    aft(SEXP list) : args(as<List>(list)) {
      init = as<vec>(args["init"]);
      X = as<mat>(args["X"]);
      XD = as<mat>(args["XD"]);
      event = as<vec>(args["event"]);
      time = as<vec>(args["time"]);
      boundaryKnots = as<vec>(args["boundaryKnots"]);
      interiorKnots = as<vec>(args["interiorKnots"]);
      q_matrix = as<mat>(args["q.const"]);
      s = ns(boundaryKnots, interiorKnots, q_matrix, 1);
      delayed = as<bool>(args["delayed"]);
      if (delayed) {
	time0 = as<vec>(args["time0"]);
	X0 = as<mat>(args["X0"]);
      }
    }
    mat rmult(mat m, vec v) {
      mat out(m);
      out.each_col() %= v;
      return out;
    }
    mat rmult(mat m, uvec v) {
      mat out(m);
      out.each_col() %= conv_to<vec>::from(v);
      return out;
    }
    double objective(vec betafull)
    {
      vec beta = betafull.subvec(0,X.n_cols-1);
      vec betas = betafull.subvec(X.n_cols,betafull.size()-1);
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec logtstar = log(time) - eta;
      vec etas = s.basis(logtstar) * betas;
      vec etaDs = s.basis(logtstar,1) * betas;
      // fix bounds on etaDs
      vec eps = etaDs*0. + 1e-8;
      double pen = dot(min(etaDs,eps), min(etaDs,eps));
      etaDs = max(etaDs, eps);
      // fix bounds on etaD
      pen += dot(min(1/time-etaD,eps), min(1/time-etaD,eps));
      etaD = 1/time - max(1/time-etaD, eps);
      vec logh = etas + log(etaDs) + log(1/time -etaD);
      vec H = exp(etas);
      double f = pen - (dot(logh,event) - sum(H));
      if (delayed) {
	vec eta0 = X0 * beta;
	vec logtstar0 = log(time0) - eta0;
	vec etas0 = s.basis(logtstar0) * betas;
	vec etaDs0 = s.basis(logtstar0,1) * betas;
	vec H0 = exp(etas0);
	vec eps0 = etaDs0*0. + 1e-8;
	pen += dot(min(etaDs0,eps0), min(etaDs0,eps0));
	f -= sum(H0);
      }
      return f;
    }
    vec gradient(vec betafull)
    {
      vec beta = betafull.subvec(0,X.n_cols-1);
      vec betas = betafull.subvec(X.n_cols,betafull.size()-1);
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec logtstar = log(time) - eta;
      mat Xs = s.basis(logtstar);
      mat XDs = s.basis(logtstar,1);
      mat XDDs = s.basis(logtstar,2);
      vec etas = Xs * betas;
      vec etaDs = XDs * betas;
      vec etaDDs = XDDs * betas;
      // H calculations
      vec H = exp(etas);
      mat dHdbetas = rmult(Xs,H);
      mat dHdbeta = -rmult(X,H % etaDs);
      // penalties
      vec eps = etaDs*0. + 1e-8;
      uvec pindexs = (etaDs < eps);
      uvec pindex = ((1.0/time - etaD) < eps);
      // fix bounds on etaDs
      mat pgrads = join_rows(-2*rmult(X,etaDs % etaDDs),2*rmult(XDs,etaDs));
      etaDs = max(etaDs, eps);
      // fix bounds on etaD
      mat pgrad = join_rows(-2*rmult(XD,1/time-etaD),XDs*0.0);
      etaD = 1/time - max(1/time-etaD, eps);
      vec logh = etas + log(etaDs) + log(1/time -etaD);
      vec h = exp(logh);
      mat dloghdbetas = Xs+rmult(XDs,1/etaDs % (1-pindexs));
      mat dloghdbeta = -rmult(X,etaDs % (1-pindexs) % (1-pindex)) - rmult(X,etaDDs/etaDs % (1-pindexs) % (1-pindex)) - rmult(XD, (1-pindexs) % (1-pindex)/(1/time-etaD));
      mat gradi = join_rows(-rmult(dloghdbeta,event)+dHdbeta, -rmult(dloghdbetas,event)+dHdbetas) + rmult(pgrad,pindex) + rmult(pgrads,pindexs);
      vec out = sum(gradi,0).t();
      if (delayed) {
	vec eta0 = X0 * beta;
	vec logtstar0 = log(time0) - eta0;
	mat Xs0 = s.basis(logtstar0);
	mat XDs0 = s.basis(logtstar0,1);
	mat XDDs0 = s.basis(logtstar0,2);
	vec etas0 = Xs0 * betas;
	vec etaDs0 = XDs0 * betas;
	vec etaDDs0 = XDDs0 * betas;
	vec H0 = exp(etas0);
	mat dHdbetas0 = rmult(Xs0,H0);
	mat dHdbeta0 = -rmult(X0,H0 % etaDs0);
	vec eps0 = etaDs0*0. + 1e-8;
	uvec pindexs0 = (etaDs0 < eps0);
	etaDs0 = max(etaDs0, eps0);
	mat pgrads0 = join_rows(-2*rmult(X0,etaDs0 % etaDDs0),2*rmult(XDs0,etaDs0));
	out += sum(join_rows(-dHdbeta0, -dHdbetas0) + rmult(pgrads0,pindexs0), 0).t();
      }
      return out;
    }
    double objective(NumericVector betafull) {
      return objective(as<vec>(wrap(betafull)));
    }
    NumericVector gradient(NumericVector betafull) {
      vec value = gradient(as<vec>(wrap(betafull)));
      return as<NumericVector>(wrap(value));
    }
    vec survival(vec time, mat X) {
      vec beta = init.subvec(0,X.n_cols-1);
      vec betas = init.subvec(X.n_cols,init.size()-1);
      vec eta = X * beta;
      vec logtstar = log(time) - eta;
      vec etas = s.basis(logtstar) * betas;
      vec S = exp(-exp(etas));
      return S;
    }
  };
  
  RcppExport SEXP aft_model_output(SEXP args) {
    aft model(args);
    List list = as<List>(args);
    std::string return_type = as<std::string>(list["return_type"]);
    if (return_type == "nmmin") {
      // model.pre_process();
      NelderMead nm;
      nm.trace = as<int>(list["trace"]);
      NumericVector betafull = as<NumericVector>(wrap(model.init));
      nm.optim<aft>(betafull,model);
      // model.post_process();
      return List::create(_("fail")=nm.fail, 
			  _("coef")=wrap(nm.coef),
			  _("hessian")=wrap(nm.hessian));
    }
    else if (return_type == "vmmin") {
      // model.pre_process();
      BFGS bfgs;
      bfgs.trace = as<int>(list["trace"]);
      NumericVector betafull = as<NumericVector>(wrap(model.init));
      bfgs.optim<aft>(betafull,model);
      // model.post_process();
      return List::create(_("fail")=bfgs.fail, 
			  _("coef")=wrap(bfgs.coef),
			  _("hessian")=wrap(bfgs.hessian));
    }
    else if (return_type == "objective")
      return wrap(model.objective(model.init));
    else if (return_type == "gradient")
      return wrap(model.gradient(model.init));
    else if (return_type == "survival")
      return wrap(model.survival(as<vec>(list["time"]),as<mat>(list["X"])));
    else {
      REprintf("Unknown return_type.\n");
      return wrap(-1);
    }
  }
  
  // RcppExport SEXP aft_objective_function(SEXP args)
  // {
  //   List list = as<List>(args);
  //   vec beta = as<vec>(list["beta"]);
  //   vec betas = as<vec>(list["betas"]);
  //   mat X = as<mat>(list["X"]);
  //   mat XD = as<mat>(list["XD"]);
  //   vec event = as<vec>(list["event"]);
  //   vec time = as<vec>(list["time"]);
  //   vec boundaryKnots = as<vec>(list["boundaryKnots"]);
  //   vec interiorKnots = as<vec>(list["interiorKnots"]);
  //   ns s(boundaryKnots, interiorKnots, 1);
  //   vec eta = X * beta;
  //   vec etaD = XD * beta;
  //   vec logtstar = log(time) - eta;
  //   vec etas = s.basis(logtstar) * betas;
  //   vec etaDs = s.basis(logtstar,1) * betas;
  //   vec eps = etaDs*0. + 1e-8;
  //   double pen = dot(min(etaDs,eps), min(etaDs,eps));
  //   etaDs = max(etaDs, eps);
  //   vec logh = etas + log(etaDs) + log(etaD+1.) - log(time);
  //   vec H = exp(etas);
  //   double f = pen - (dot(logh,event) - sum(H));
  //   return wrap(f);
  // }

} // namespace rstpm2

