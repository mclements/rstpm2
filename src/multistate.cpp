#include <RcppArmadillo.h>

using namespace arma;
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
// Rcpp::List multistate_ddt(mat P, cube Pu, cube Q, cube Qu) {
RcppExport SEXP multistate_ddt(SEXP P_, SEXP Pu_, SEXP Q_, SEXP Qu_) {
  // P(nstates,nobs), Pu(nstates,ncoef,nobs), Q(nstates,nstates,nobs), Qu(nstates*nstates,ncoef,nobs)
  mat P = Rcpp::as<mat>(P_);
  cube Pu = Rcpp::as<cube>(Pu_);
  cube Q = Rcpp::as<cube>(Q_);
  cube Qu = Rcpp::as<cube>(Qu_);
  int nstates = P.n_rows;
  int nobs=P.n_cols;
  int ncoef=Pu.n_cols;
  mat dPdt = P*0.0;   // nstates,nobs
  cube dPudt = Pu*0.0; // nstates,ncoef,nobs
  for(int i=0; i<nobs; i++) {
    mat Qi = Q.slice(i);
    rowvec Pi = conv_to<rowvec>::from(P.col(i));
    vec rs = sum(Qi,1);
    for (int j=0; j<nstates; j++)
      Qi(j,j) = -rs(j);
    dPdt.col(i) = conv_to<vec>::from(Pi*Qi); // is this valid?
    cube Qui = Qu(span::all, span::all, span(i,i));
    Qui.reshape(nstates,nstates,ncoef);
    mat Pui = Pu.slice(i);
    for (int j=0; j<ncoef; j++) {
      vec rs = sum(Qui.slice(j),1);
      for (int m=0; m<nstates; m++)
	Qui(m,m,j) = -rs(m);
    }
    dPudt.slice(i) = (Pui.t() * Qi).t();
    for(int j=0; j<ncoef; j++) {
      rowvec PiQuij = Pi * Qui.slice(j);
      for(int k=0; k<nstates; k++)
	dPudt(k,j,i) += PiQuij(k);
    }
  }
  return Rcpp::wrap(Rcpp::List::create(dPdt,dPudt));
}
