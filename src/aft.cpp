#include <RcppArmadillo.h>
#include <c_optim.h>
#include <splines.h>

namespace rstpm2 {
  using namespace arma;
  using namespace Rcpp;

  class aft {
  public:
    List args;
    vec init;
    mat X, X0;
    mat XD, XD0;
    vec event;
    vec time, time0;
    vec boundaryKnots;
    vec interiorKnots;
    mat q_matrix;
    int cure;
    ns s;
    bool delayed;
    double kappa; // scale for the quadratic penalty for monotone splines
    aft(SEXP list) : args(as<List>(list)) {
      init = as<vec>(args["init"]);
      X = as<mat>(args["X"]);
      XD = as<mat>(args["XD"]);
      XD0 = as<mat>(args["XD0"]);
      event = as<vec>(args["event"]);
      time = as<vec>(args["time"]);
      boundaryKnots = as<vec>(args["boundaryKnots"]);
      interiorKnots = as<vec>(args["interiorKnots"]);
      q_matrix = as<mat>(args["q.const"]);
      cure = as<int>(args["cure"]);
      s = ns(boundaryKnots, interiorKnots, q_matrix, 1, cure);
      delayed = as<bool>(args["delayed"]);
      if (delayed) {
	time0 = as<vec>(args["time0"]);
	X0 = as<mat>(args["X0"]);
      }
      kappa = 1.0e3;
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
      // add penalty for monotone splines
      vec betasStar = s.q_matrix.t() * betas;
      for (size_t i=1; i<betasStar.size(); i++) {
      	double delta = betasStar(i)-betasStar(i-1);
      	if (delta<0.0)
      	  pen += kappa*delta*delta;
      }
      vec logh = etas + log(etaDs) + log(1/time -etaD);
      vec H = exp(etas);
      double f = pen - (dot(logh,event) - sum(H));
      if (delayed) {
	vec eta0 = X0 * beta;
	vec logtstar0 = log(time0) - eta0;
	vec etas0 = s.basis(logtstar0) * betas;
	vec etaDs0 = s.basis(logtstar0,1) * betas;
	vec H0 = exp(etas0);
	vec eps0 = etaDs0*0. + 1e-16;
	f += dot(min(etaDs0,eps0), min(etaDs0,eps0));
	f -= sum(H0);
      }
      return f;
    }
    vec gradientPenalty(mat Q, vec beta) { // Q: (nbeta+2) x nbeta
      size_t n = Q.n_rows;
      mat D = join_rows(zeros(n-1,1),eye(n-1,n-1)) - join_rows(eye(n-1,n-1),zeros(n-1,1)); // (nbeta+1) x (nbeta+2)
      vec delta = D * Q * beta; // nbeta+1
      mat M = Q.t() * D.row(0).t() * D.row(0) * Q * (delta(0)<0.0); // nbeta x nbeta
      for(size_t j=1; j<delta.size(); j++) {
	if (delta(j)<0.0)
	  M += Q.t() * D.row(j).t() * D.row(j) * Q;
      }
      return 2*M*beta*kappa;
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
      out += join_cols(beta*0.0, gradientPenalty(s.q_matrix.t(), betas));
      if (delayed) {
	vec eta0 = X0 * beta;
	vec etaD0 = XD0 * beta;
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
	uvec pindex0 = ((1.0/time0 - etaD0) < eps0);
	uvec pindexs0 = (etaDs0 < eps0);
	etaDs0 = max(etaDs0, eps0);
	mat pgrad0 = join_rows(-2*rmult(XD0,1/time0-etaD0),XDs0*0.0);
	mat pgrads0 = join_rows(-2*rmult(X0,etaDs0 % etaDDs0),2*rmult(XDs0,etaDs0));
	out += sum(join_rows(-dHdbeta0, -dHdbetas0) + rmult(pgrads0,pindexs0) +
		   rmult(pgrad0,pindex0), 0).t();
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

    mat haz(vec time, mat X, mat XD)
    {
      vec beta = init.subvec(0,X.n_cols-1);
      vec betas = init.subvec(X.n_cols,init.size()-1);
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec logtstar = log(time) - eta;
      mat Xs = s.basis(logtstar);
      mat XDs = s.basis(logtstar,1);
      mat XDDs = s.basis(logtstar,2);
      vec etas = Xs * betas;
      vec etaDs = XDs * betas;
      vec etaDDs = XDDs * betas;
      // penalties
      vec eps = etaDs*0. + 1e-8;
      uvec pindexs = (etaDs < eps);
      uvec pindex = ((1.0/time - etaD) < eps);
      // fix bounds on etaDs
      etaDs = max(etaDs, eps);
      // fix bounds on etaD
      etaD = 1/time - max(1/time-etaD, eps);
      vec logh = etas + log(etaDs) + log(1/time -etaD);
      vec h = exp(logh);
      return h;
    }
  
    mat gradh(vec time, mat X, mat XD)
    {
      vec beta = init.subvec(0,X.n_cols-1);
      vec betas = init.subvec(X.n_cols,init.size()-1);
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec logtstar = log(time) - eta;
      mat Xs = s.basis(logtstar);
      mat XDs = s.basis(logtstar,1);
      mat XDDs = s.basis(logtstar,2);
      vec etas = Xs * betas;
      vec etaDs = XDs * betas;
      vec etaDDs = XDDs * betas;
      // penalties
      vec eps = etaDs*0. + 1e-8;
      uvec pindexs = (etaDs < eps);
      uvec pindex = ((1.0/time - etaD) < eps);
      // fix bounds on etaDs
      etaDs = max(etaDs, eps);
      // fix bounds on etaD
      etaD = 1/time - max(1/time-etaD, eps);
      vec logh = etas + log(etaDs) + log(1/time -etaD);
      vec h = exp(logh);
      mat dloghdbetas = Xs+rmult(XDs,1/etaDs % (1-pindexs));
      mat dloghdbeta = -rmult(X,etaDs % (1-pindexs) % (1-pindex)) - rmult(X,etaDDs/etaDs % (1-pindexs) % (1-pindex)) - rmult(XD, (1-pindexs) % (1-pindex)/(1/time-etaD));
      mat gradh = join_rows(rmult(dloghdbeta,h), rmult(dloghdbetas,h));
      return gradh;
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
      nm.maxit = as<int>(list["maxit"]);
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
      bfgs.maxit = as<int>(list["maxit"]);
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
    else if (return_type == "haz")
      return wrap(model.haz(as<vec>(list["time"]),as<mat>(list["X"]),as<mat>(list["XD"])));
    else if (return_type == "gradh")
      return wrap(model.gradh(as<vec>(list["time"]),as<mat>(list["X"]),as<mat>(list["XD"])));
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

