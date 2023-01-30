#include <RcppArmadillo.h>
#include <c_optim.h>
#include <splines.h>
#include <vector>

namespace rstpm2 {
  using namespace arma;
  using namespace Rcpp;


  // TODO:
  // include mixture models in .gradh(), .haz() and .survival()
  
  class aft_integrated {
  public:
    List args;
    vec init;
    std::vector<mat> X_vector, X_vector0;
    mat Xc;
    mat Xt;
    vec event;
    vec time, time0;
    vec gnodes, gweights;
    vec boundaryKnots;
    vec interiorKnots;
    mat q_matrix;
    int cure;
    int mixture;
    ns s;
    double kappa; // scale for the quadratic penalty for monotone splines
    uvec delayed;
    aft_integrated(Rcpp::List args) : args(args) {
      init = as<vec>(args("init"));
      X_vector = as<std::vector<mat>>(args("X_list"));
      X_vector0 = as<std::vector<mat>>(args("X_list0"));
      Xt = as<mat>(args("Xt"));
      Xc = as<mat>(args("Xc"));
      gnodes = as<vec>(args("gnodes"));
      gweights = as<vec>(args("gweights"));
      event = as<vec>(args("event"));
      time = as<vec>(args("time"));
      time0 = as<vec>(args("time0"));
      boundaryKnots = as<vec>(args("boundaryKnots"));
      interiorKnots = as<vec>(args("interiorKnots"));
      q_matrix = as<mat>(args("q.const"));
      cure = as<int>(args("cure"));
      mixture = as<int>(args("mixture"));
      s = ns(boundaryKnots, interiorKnots, q_matrix, 1, cure);
      kappa = 1.0e3;
      delayed = time0>0.0;
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
	beta = betafull.subvec(0,Xt.n_cols-1);
	betac = betafull.subvec(Xt.n_cols, Xt.n_cols + Xc.n_cols - 1);
	betas = betafull.subvec(Xt.n_cols+Xc.n_cols, betafull.size()-1);
      } else {
	beta = betafull.subvec(0,Xt.n_cols-1);
	betas = betafull.subvec(Xt.n_cols, betafull.size()-1);
      }
      vec tstar = time*0.0;
      vec scale = time/2.0;
      // Can this be done more efficiently using a cube?
      for(size_t i=0; i<X_vector.size(); i++) {
	tstar += gweights[i]*scale % exp(-X_vector[i]*beta);
      }
      vec logtstar = log(tstar);
      vec etat = Xt*beta;
      vec etas = s.basis(logtstar) * betas;
      vec etaDs = s.basis(logtstar,1) * betas;
      // fix bounds on etaDs
      vec eps = etaDs*0. + 1e-8;
      double pen = dot(min(etaDs,eps), min(etaDs,eps));
      etaDs = max(etaDs, eps);
      // add penalty for monotone splines
      vec betasStar = s.q_matrix.t() * betas;
      for (size_t i=1; i<betasStar.size(); i++) {
      	double delta = betasStar(i)-betasStar(i-1);
      	if (delta<0.0)
      	  pen += kappa*delta*delta;
      }
      vec logh = etas + log(etaDs) - etat - logtstar;
      vec H = exp(etas);
      double ll = dot(logh,event) - sum(H) - pen; // NB: addition in the gradient
      if (mixture) {
	etac = Xc * betac;
        vec cure_frac = exp(etac)/(1.0+exp(etac));
	vec Hu = H, loghu=logh;
	vec S = cure_frac + (1-cure_frac) % exp(-Hu);
        ll = dot(event, log(1-cure_frac)+loghu-Hu) +
	  dot(1-event,log(S)) - pen;
      }
      if (any(delayed)) {
	vec tstar0 = time0(delayed)*0.0;
	vec scale0 = time0(delayed)/2.0;
	for(size_t i=0; i<X_vector0.size(); i++) {
	  mat X0 = X_vector0[i].rows(delayed);
	  tstar0 += gweights[i]*scale0 % exp(-X0*beta);
	}
	vec logtstar0 = log(tstar0);
	vec etas0 = s.basis(logtstar0) * betas;
	vec etaDs0 = s.basis(logtstar0,1) * betas;
	vec H0 = exp(etas0);
	vec eps0 = etaDs0*0. + 1e-16;
	ll -= dot(min(etaDs0,eps0), min(etaDs0,eps0));
	if (!mixture) {
	  ll += sum(H0);
	} else {
	  vec cure_frac0 = exp(etac(delayed))/(1+exp(etac(delayed)));
	  vec S0_mix = cure_frac0 + (1.0-cure_frac0) % exp(-H0);
	  vec H0_mix = -log(S0_mix);
	  ll += sum(H0_mix);
	}
	// if (pen>0.0) ll=-100;
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
	beta = betafull.subvec(0,Xt.n_cols-1);
	betac = betafull.subvec(Xt.n_cols, Xt.n_cols + Xc.n_cols - 1);
	betas = betafull.subvec(Xt.n_cols+Xc.n_cols, betafull.size()-1);
      } else {
	beta = betafull.subvec(0,Xt.n_cols-1);
	betas = betafull.subvec(Xt.n_cols, betafull.size()-1);
      }
      vec tstar = time*0.0;
      mat Xtstar = Xt*0.0;
      vec scale = time/2.0;
      // Can this be done more efficiently using a cube?
      for(size_t i=0; i<X_vector.size(); i++) {
	vec integrand = gweights[i]*scale % exp(-X_vector[i]*beta);
	tstar += integrand;
	Xtstar += rmult(X_vector[i],integrand);
      }
      vec logtstar = log(tstar);
      vec etat = Xt*beta;
      mat Xs = s.basis(logtstar);
      mat XDs = s.basis(logtstar,1);
      mat XDDs = s.basis(logtstar,2);
      vec etas = Xs * betas;
      vec etaDs = XDs * betas;
      vec etaDs_old = etaDs;
      vec etaDDs = XDDs * betas;
      // H calculations
      vec H = exp(etas);
      mat dHdbeta = -rmult(Xtstar, H % etaDs / tstar);
      // penalties
      vec eps = etaDs*0. + 1e-8;
      uvec pindexs = (etaDs < eps);
      // fix bounds on etaDs
      mat pgrad = join_rows(2*rmult(Xtstar, etaDs % etaDDs / tstar), -2.0*rmult(XDs,etaDs));
      etaDs = max(etaDs, eps);
      mat dHdbetas = rmult(Xs,H);
      vec logh = etas + log(etaDs) - etat - logtstar;
      vec h = exp(logh);
      mat dloghdbetas = Xs+rmult(XDs,1/etaDs % (1-pindexs));
      mat dloghdbeta = - rmult(Xtstar,etaDDs/etaDs_old/tstar) -
	rmult(Xtstar,etaDs_old / tstar) + rmult(Xtstar,1.0/tstar) - Xt;
      mat gradi = join_rows(rmult(dloghdbeta,event)-dHdbeta,
			    rmult(dloghdbetas,event)-dHdbetas) + rmult(pgrad,pindexs);
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
	mat pgrad = join_rows(2*rmult(Xtstar, etaDs_old % etaDDs / tstar), Xc*0.0, -2.0*rmult(XDs,etaDs_old));
	gradi = join_rows(rmult(dloghdbeta_mix,event)-dHdbeta_mix,
			  rmult(dloghdtheta_mix,event)-dHdtheta_mix,
			  rmult(dloghdbetas_mix,event)-dHdbetas_mix) + rmult(pgrad,pindexs);
	gr = sum(gradi,0).t();
	gr -= join_cols(beta*0.0, betac*0.0, gradientPenalty(s.q_matrix.t(), betas));
      }
      if (any(delayed)) {
	vec tstar0 = time0(delayed)*0.0;
	vec scale0 = time0(delayed)/2.0;
	mat Xtstar0 = Xt.rows(delayed)*0.0;
	// Can this be done more efficiently using a cube?
	for(size_t i=0; i<X_vector0.size(); i++) {
	  vec integrand = gweights[i]*scale0 % exp(-X_vector0[i].rows(delayed)*beta);
	  tstar0 += integrand;
	  Xtstar0 += rmult(X_vector0[i].rows(delayed),integrand);
	}
	vec logtstar0 = log(tstar0);
	mat Xs0 = s.basis(logtstar0);
	mat XDs0 = s.basis(logtstar0,1);
	mat XDDs0 = s.basis(logtstar0,2);
	vec etas0 = Xs0 * betas;
	vec etaDs0 = XDs0 * betas;
	vec etaDs0_old = etaDs0;
	vec etaDDs0 = XDDs0 * betas;
	vec H0 = exp(etas0);
	mat dHdbetas0 = rmult(Xs0,H0);
	mat dHdbeta0 = -rmult(Xtstar0,H0 % etaDs0 / tstar0);
	vec eps0 = etaDs0*0. + 1e-8;
	uvec pindexs0 = (etaDs0 < eps0);
	etaDs0 = max(etaDs0, eps0);
	mat pgrad0 = join_rows(2*rmult(Xtstar0, etaDs0_old % etaDDs0 / tstar0), Xs*0.0);
	mat pgrads0 = join_rows(Xt*0.0, -2.0*rmult(XDs0,etaDs0_old));
	if (!mixture) {
	  gr += sum(join_rows(dHdbeta0, dHdbetas0) + rmult(pgrads0,pindexs0) + rmult(pgrad0,pindexs0), 0).t();
	} else {
	  vec cure_frac0 = exp(etac(delayed))/(1.0+exp(etac(delayed)));
	  vec S_mix0 = cure_frac0 + (1.0-cure_frac0) % exp(-H0);
	  vec H_mix0 = -log(S_mix0);
	  mat dpidtheta0 = rmult(Xc.rows(delayed), cure_frac0 % (1-cure_frac0));
	  mat dHdbeta_mix0 = rmult(dHdbeta0,(1.0 - cure_frac0) % exp(-H0) / S_mix0);
	  mat dHdbetas_mix0 = rmult(dHdbetas0,(1.0 - cure_frac0) % exp(-H0) / S_mix0);
	  mat dHdtheta_mix0 = -rmult(dpidtheta0, (1.0 - exp(-H0))/S_mix0);
	  mat pgrad0 = join_rows(2*rmult(Xtstar0, etaDs0_old % etaDDs0 / tstar0), Xc*0.0, Xs*0.0);
	  mat pgrads0 = join_rows(Xt*0.0, Xc*0.0, -2.0*rmult(XDs0,etaDs0_old));
	  gr += sum(join_rows(dHdbeta_mix0, dHdtheta_mix0, dHdbetas_mix0) +
		    rmult(pgrads0,pindexs0) + rmult(pgrad0,pindexs0), 0).t();
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
    // vec survival(vec time, std::vector<mat> X_vector1) {
    //   vec beta = init.subvec(0,X.n_cols-1);
    //   vec betas = init.subvec(X.n_cols,init.size()-1);
    //   vec eta = X * beta;
    //   vec logtstar = log(time) - eta;
    //   vec etas = s.basis(logtstar) * betas;
    //   vec S = exp(-exp(etas));
    //   return S;
    // }
    // mat haz(vec time, mat X, mat XD)
    // {
    //   vec beta = init.subvec(0,X.n_cols-1);
    //   vec betas = init.subvec(X.n_cols,init.size()-1);
    //   vec eta = X * beta;
    //   vec etaD = XD * beta;
    //   vec logtstar = log(time) - eta;
    //   mat Xs = s.basis(logtstar);
    //   mat XDs = s.basis(logtstar,1);
    //   mat XDDs = s.basis(logtstar,2);
    //   vec etas = Xs * betas;
    //   vec etaDs = XDs * betas;
    //   vec etaDDs = XDDs * betas;
    //   // penalties
    //   vec eps = etaDs*0. + 1e-8;
    //   uvec pindexs = (etaDs < eps);
    //   uvec pindex = ((1.0/time - etaD) < eps);
    //   // fix bounds on etaDs
    //   etaDs = max(etaDs, eps);
    //   // fix bounds on etaD
    //   etaD = 1/time - max(1/time-etaD, eps);
    //   vec logh = etas + log(etaDs) + log(1/time -etaD);
    //   vec h = exp(logh);
    //   return h;
    // }
    // mat gradh(vec time, mat X, mat XD)
    // {
    //   vec beta = init.subvec(0,X.n_cols-1);
    //   vec betas = init.subvec(X.n_cols,init.size()-1);
    //   vec eta = X * beta;
    //   vec etaD = XD * beta;
    //   vec logtstar = log(time) - eta;
    //   mat Xs = s.basis(logtstar);
    //   mat XDs = s.basis(logtstar,1);
    //   mat XDDs = s.basis(logtstar,2);
    //   vec etas = Xs * betas;
    //   vec etaDs = XDs * betas;
    //   vec etaDDs = XDDs * betas;
    //   // penalties
    //   vec eps = etaDs*0. + 1e-8;
    //   uvec pindexs = (etaDs < eps);
    //   uvec pindex = ((1.0/time - etaD) < eps);
    //   // fix bounds on etaDs
    //   etaDs = max(etaDs, eps);
    //   // fix bounds on etaD
    //   etaD = 1/time - max(1/time-etaD, eps);
    //   vec logh = etas + log(etaDs) + log(1/time -etaD);
    //   vec h = exp(logh);
    //   mat dloghdbetas = Xs+rmult(XDs,1/etaDs % (1-pindexs));
    //   mat dloghdbeta = -rmult(X,etaDs % (1-pindexs) % (1-pindex)) - rmult(X,etaDDs/etaDs % (1-pindexs) % (1-pindex)) - rmult(XD, (1-pindexs) % (1-pindex)/(1/time-etaD));
    //   mat gradh = join_rows(rmult(dloghdbeta,h), rmult(dloghdbetas,h));
    //   return gradh;
    // }
  };
  

    RcppExport SEXP aft_integrated_model_output(SEXP args) {
    using namespace Rcpp;
    using namespace arma;
    aft_integrated model(args);
    List list = as<List>(args);
    std::string return_type = as<std::string>(list["return_type"]);
    if (return_type == "nmmin") {
      // model.pre_process();
      NelderMead nm;
      nm.trace = as<int>(list["trace"]);
      nm.maxit = as<int>(list["maxit"]);
      NumericVector betafull = as<NumericVector>(wrap(model.init));
      nm.optim<aft_integrated>(betafull,model);
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
      bfgs.optim<aft_integrated>(betafull,model);
      // model.post_process();
      return List::create(_("fail")=bfgs.fail, 
			  _("coef")=wrap(bfgs.coef),
			  _("hessian")=wrap(bfgs.hessian));
    }
    else if (return_type == "objective")
      return wrap(model.objective(model.init));
    else if (return_type == "gradient")
      return wrap(model.gradient(model.init));
    // else if (return_type == "survival")
    //   return wrap(model.survival(as<vec>(list["time"]),as<mat>(list["X"])));
    // else if (return_type == "haz")
    //   return wrap(model.haz(as<vec>(list["time"]),as<mat>(list["X"]),as<mat>(list["XD"])));
    // else if (return_type == "gradh")
    //   return wrap(model.gradh(as<vec>(list["time"]),as<mat>(list["X"]),as<mat>(list["XD"])));
    else {
      REprintf("Unknown return_type.\n");
      return wrap(-1);
    }
  }

} // namespace rstpm2


