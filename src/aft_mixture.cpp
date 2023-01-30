#include <RcppArmadillo.h>
#include <c_optim.h>
#include <splines.h>
#include <vector>

namespace rstpm2 {
  using namespace arma;
  using namespace Rcpp;


  // TODO:
  // include mixture models in .gradh(), .haz() and .survival()
  
  class aft_mixture {
  public:
    List args;
    vec init;
    mat X, X0, Xc, Xc0;
    mat XD, XD0;
    vec event;
    vec time, time0;
    vec boundaryKnots;
    vec interiorKnots;
    mat q_matrix;
    int cure;
    int mixture;
    ns s;
    bool delayed;
    double kappa; // scale for the quadratic penalty for monotone splines
    aft_mixture(Rcpp::List args) : args(args) {
      init = as<vec>(args("init"));
      X = as<mat>(args("X"));
      Xc = as<mat>(args("Xc"));
      Xc0 = as<mat>(args("Xc0"));
      XD = as<mat>(args("XD"));
      XD0 = as<mat>(args("XD0"));
      event = as<vec>(args("event"));
      time = as<vec>(args("time"));
      boundaryKnots = as<vec>(args("boundaryKnots"));
      interiorKnots = as<vec>(args("interiorKnots"));
      q_matrix = as<mat>(args("q.const"));
      cure = as<int>(args("cure"));
      mixture = as<int>(args("mixture"));
      s = ns(boundaryKnots, interiorKnots, q_matrix, 1, cure);
      delayed = as<bool>(args("delayed"));
      if (delayed) {
	time0 = as<vec>(args("time0"));
	X0 = as<mat>(args("X0"));
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
      vec beta, betac, betas, etac;
      if (mixture) {
	beta = betafull.subvec(0,X.n_cols-1);
	betac = betafull.subvec(X.n_cols, X.n_cols + Xc.n_cols - 1);
	betas = betafull.subvec(X.n_cols+Xc.n_cols, betafull.size()-1);
      } else {
	beta = betafull.subvec(0,X.n_cols-1);
	betas = betafull.subvec(X.n_cols, betafull.size()-1);
      }
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec logtstar = log(time) - eta;
      vec etas = s.basis(logtstar) * betas;
      vec etaDs = s.basis(logtstar,1) * betas;
      vec etaDs_old = etaDs;
      // fix bounds on etaDs
      vec eps = etaDs*0. + 1e-8;
      double pen = dot(min(etaDs,eps), min(etaDs,eps));
      etaDs = max(etaDs, eps);
      // fix bounds on etaD
      pen += dot(min(1-etaD,eps), min(1-etaD,eps));
      etaD = 1 - max(1-etaD, eps);
      // add penalty for monotone splines
      vec betasStar = s.q_matrix.t() * betas;
      for (size_t i=1; i<betasStar.size(); i++) {
      	double delta = betasStar(i)-betasStar(i-1);
      	if (delta<0.0)
      	  pen += kappa*delta*delta;
      }
      vec logh = etas + log(etaDs) + log(1/time -etaD/time);
      vec H = exp(etas);
      double ll = dot(logh,event) - sum(H) - pen;
      if (mixture) {
	etac = Xc * betac;
        vec cure_frac = exp(etac)/(1+exp(etac));
        ll = dot(event, log(1-cure_frac)+logh-H)
	  + dot(1-event,log(cure_frac+(1-cure_frac)%exp(-H))) - pen;
      }
      if (delayed) {
	vec eta0 = X0 * beta;
	vec logtstar0 = log(time0) - eta0;
	vec etas0 = s.basis(logtstar0) * betas;
	vec etaDs0 = s.basis(logtstar0,1) * betas;
	vec etaDs0_old = etaDs0;
	vec H0 = exp(etas0);
	vec eps0 = etaDs0*0. + 1e-16;
	ll -= dot(min(etaDs0,eps0), min(etaDs0,eps0));
	if (!mixture) {
	  ll += sum(H0);
	} else {
	  vec etac0 = Xc0 * betac;
	  vec cure_frac0 = exp(etac0)/(1+exp(etac0));
	  vec S0_mix = cure_frac0 + (1.0-cure_frac0) % exp(-H0);
	  vec H0_mix = -log(S0_mix);
	  ll += sum(H0_mix);
	}
      }
      return -ll;
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
      vec beta, betac, betas, etac;
      if (mixture) {
	beta = betafull.subvec(0,X.n_cols-1);
	betac = betafull.subvec(X.n_cols, X.n_cols + Xc.n_cols - 1);
	betas = betafull.subvec(X.n_cols+Xc.n_cols, betafull.size()-1);
      } else {
	beta = betafull.subvec(0,X.n_cols-1);
	betas = betafull.subvec(X.n_cols, betafull.size()-1);
      }
      vec eta = X * beta;
      vec etaD = XD * beta;
      vec etaD_old = etaD;
      vec logtstar = log(time) - eta;
      vec tstar = exp(logtstar);
      mat Xs = s.basis(logtstar);
      mat XDs = s.basis(logtstar,1);
      mat XDDs = s.basis(logtstar,2);
      vec etas = Xs * betas;
      vec etaDs = XDs * betas;
      vec etaDs_old = etaDs;
      vec etaDDs = XDDs * betas;
      // H calculations
      vec H = exp(etas);
      mat dHdbeta = -rmult(X,H % etaDs);
      mat dHdbetas = rmult(Xs,H);
      // penalties
      vec eps = etaDs*0. + 1e-8;
      uvec pindexs = (etaDs < eps);
      uvec pindex = ((1.0 - etaD) < eps);
      // fix bounds on etaDs
      mat pgrads = join_rows(X*0.0,-2*rmult(XDs,etaDs));
      // fix bounds on etaD
      mat pgrad = join_rows(2.0*rmult(X,etaDs % etaDDs)+2.0*rmult(XD,1-etaD),Xs*0.0);
      etaDs = max(etaDs, eps);
      etaD = 1 - max(1-etaD, eps);
      vec logh = etas + log(etaDs) + log(1/time -etaD/time);
      vec h = exp(logh);
      mat dloghdbetas = Xs+rmult(XDs,1/etaDs % (1-pindexs));
      mat dloghdbeta = -rmult(X,etaDDs/etaDs % (1-pindexs)) - rmult(X,etaDs_old) - rmult(XD, 1/(1-etaD_old) % (1-pindex));
      mat gradi = join_rows(rmult(dloghdbeta,event)-dHdbeta, rmult(dloghdbetas,event)-dHdbetas) +
	rmult(pgrad,pindex) + rmult(pgrads,pindexs);
      vec gr = sum(gradi,0).t();
      gr -= join_cols(beta*0.0, gradientPenalty(s.q_matrix.t(), betas));
      if (mixture) {
	etac = Xc * betac;
        vec cure_frac = exp(etac)/(1.0+exp(etac));
	vec S_mix = cure_frac + (1.0-cure_frac) % exp(-H);
	vec H_mix = -log(S_mix);
	mat dpidtheta = rmult(Xc, cure_frac % (1-cure_frac));
	mat dHdbeta_mix = rmult(dHdbeta,(1.0 - cure_frac) % exp(-H) / S_mix);
	mat dHdbetas_mix = rmult(dHdbetas,(1.0 - cure_frac) % exp(-H) / S_mix);
	mat dHdtheta_mix = -rmult(dpidtheta, (1.0 - exp(-H))/S_mix);
	mat dloghdbeta_mix = dloghdbeta - dHdbeta + dHdbeta_mix;
	mat dloghdbetas_mix = dloghdbetas - dHdbetas + dHdbetas_mix;
	mat dloghdtheta_mix = dHdtheta_mix - rmult(dpidtheta,1.0 / (1.0 - cure_frac));
	mat pgrads = join_rows(X*0.0,Xc*0.0,-2*rmult(XDs,etaDs_old));
	mat pgrad = join_rows(2.0*rmult(X,etaDs_old % etaDDs)+2.0*rmult(XD,1-etaD_old),Xc*0.0,Xs*0.0);
	gradi = join_rows(rmult(dloghdbeta_mix,event)-dHdbeta_mix,
			  rmult(dloghdtheta_mix,event)-dHdtheta_mix,
			  rmult(dloghdbetas_mix,event)-dHdbetas_mix) +
	  rmult(pgrad,pindex) + rmult(pgrads,pindexs);
	gr = sum(gradi,0).t();
	gr -= join_cols(beta*0.0, betac*0.0, gradientPenalty(s.q_matrix.t(), betas));
      }
      if (delayed) {
	vec eta0 = X0 * beta;
	vec etaD0 = XD0 * beta;
	vec etaD0_old = etaD0;
	vec logtstar0 = log(time0) - eta0;
	mat Xs0 = s.basis(logtstar0);
	mat XDs0 = s.basis(logtstar0,1);
	mat XDDs0 = s.basis(logtstar0,2);
	vec etas0 = Xs0 * betas;
	vec etaDs0 = XDs0 * betas;
	vec etaDs0_old = etaDs0;
	vec etaDDs0 = XDDs0 * betas;
	vec H0 = exp(etas0);
	mat dHdbetas0 = rmult(Xs0,H0);
	mat dHdbeta0 = -rmult(X0,H0 % etaDs0);
	vec eps0 = etaDs0*0. + 1e-8;
	uvec pindex0 = ((1.0/time0 - etaD0) < eps0);
	uvec pindexs0 = (etaDs0 < eps0);
	etaDs0 = max(etaDs0, eps0);
	etaD0 = 1 - max(1-etaD0, eps0);
	if (!mixture) {
	  mat pgrad0 = join_rows(2.0*rmult(X0,etaDs0_old % etaDDs0)+2.0*rmult(XD0,1-etaD0_old),Xs0*0.0);
	  mat pgrads0 = join_rows(X0*0.0,-2*rmult(XDs0,etaDs0_old));
	  gr += sum(join_rows(dHdbeta0, dHdbetas0) + rmult(pgrads0,pindexs0) +
		    rmult(pgrad0,pindex0), 0).t();
	} else {
	  vec etac0 = Xc0 * betac;
	  vec cure_frac0 = exp(etac0)/(1.0+exp(etac0));
	  vec S_mix0 = cure_frac0 + (1.0-cure_frac0) % exp(-H0);
	  vec H_mix0 = -log(S_mix0);
	  mat dpidtheta0 = rmult(Xc0, cure_frac0 % (1-cure_frac0));
	  mat dHdbeta_mix0 = rmult(dHdbeta0,(1.0 - cure_frac0) % exp(-H0) / S_mix0);
	  mat dHdbetas_mix0 = rmult(dHdbetas0,(1.0 - cure_frac0) % exp(-H0) / S_mix0);
	  mat dHdtheta_mix0 = -rmult(dpidtheta0, (1.0 - exp(-H0))/S_mix0);
	  mat pgrad0 = join_rows(2.0*rmult(X0,etaDs0_old % etaDDs0)+2.0*rmult(XD0,1-etaD0),
				 Xc0*0.0,Xs0*0.0);
	  mat pgrads0 = join_rows(X0*0.0,Xc0*0.0,-2*rmult(XDs0,etaDs0_old));
	  gr += sum(join_rows(dHdbeta_mix0, dHdtheta_mix0, dHdbetas_mix0) +
		     rmult(pgrads0,pindexs0) +
		     rmult(pgrad0,pindex0), 0).t();
	}
      }
      return -gr;
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

  RcppExport SEXP aft_mixture_model_output(SEXP args) {
    using namespace Rcpp;
    using namespace arma;
    aft_mixture model(args);
    List list = as<List>(args);
    std::string return_type = as<std::string>(list["return_type"]);
    if (return_type == "nmmin") {
      // model.pre_process();
      NelderMead nm;
      nm.trace = as<int>(list["trace"]);
      nm.maxit = as<int>(list["maxit"]);
      NumericVector betafull = as<NumericVector>(wrap(model.init));
      nm.optim<aft_mixture>(betafull,model);
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
      bfgs.optim<aft_mixture>(betafull,model);
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

} // namespace rstpm2


