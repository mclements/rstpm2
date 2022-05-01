
#include <RcppArmadillo.h>
#include <algorithm>
#include <RcppCommon.h>

using namespace arma;

/**
   Value to pass back from pluginEstimateDiscrete
 */
struct PluginEstimateDiscrete {
  mat X;           /**< State matrix (nState x nTimes) */
  mat variance;    /**< State variance matrix (nState x nTimes) */
  cube covariance; /**< State covariance matrix (nState x nState x nTimes) */
  bool vcov;       /**< Indicator whether covariances are recorded */
  int n;           /**< Number of initial observations */
  mat Y;           /**< Matrix of weighted X */
  mat varY;        /**< Matrix of variances for weighted X */
};
/**
   Value to pass back from pluginEstimateCts
 */
struct PluginEstimateCts {
  mat X;           /**< State matrix (nState x nTimes) */
  mat variance;    /**< State variance matrix (nState x nTimes) */
  cube covariance; /**< State covariance matrix (nState x nState x nTimes) */
  bool vcov;       /**< Indicator whether covariances are recorded */
  int n;           /**< Number of initial observations */
  vec time;       /**< Times */
  mat Y;           /**< Matrix of weighted X */
  mat varY;        /**< Matrix of variances for weighted X */
};

namespace Rcpp {
  template <>
  SEXP wrap(const PluginEstimateDiscrete& x);
  template <>
  SEXP wrap(const PluginEstimateCts& x);
}

#include <Rcpp.h>

namespace Rcpp {
  template <>
  SEXP wrap(const PluginEstimateDiscrete& x) {
    return Rcpp::wrap(Rcpp::List::create(Rcpp::Named("X") = Rcpp::wrap(x.X),
					 Rcpp::Named("variance") = Rcpp::wrap(x.variance),
					 Rcpp::Named("covariance") = Rcpp::wrap(x.covariance),
					 Rcpp::Named("vcov") = Rcpp::wrap(x.vcov),
					 Rcpp::Named("n") = Rcpp::wrap(x.n),
					 Rcpp::Named("Y") = Rcpp::wrap(x.Y),
					 Rcpp::Named("varY") = Rcpp::wrap(x.varY)));
  };
  template <>
  SEXP wrap(const PluginEstimateCts& x) {
    return Rcpp::wrap(Rcpp::List::create(Rcpp::Named("X") = Rcpp::wrap(x.X),
					 Rcpp::Named("variance") = Rcpp::wrap(x.variance),
					 Rcpp::Named("covariance") = Rcpp::wrap(x.covariance),
					 Rcpp::Named("vcov") = Rcpp::wrap(x.vcov),
					 Rcpp::Named("n") = Rcpp::wrap(x.n),
					 Rcpp::Named("time") = Rcpp::wrap(x.time),
					 Rcpp::Named("Y") = Rcpp::wrap(x.Y),
					 Rcpp::Named("varY") = Rcpp::wrap(x.varY)));
  };
}


/**
   Find values within an interval
 */
class FindInterval {
public:
  typedef std::vector<double> stdvec;
  typedef stdvec::iterator Iter; /**< Iterator used for speed cf. convenience */
  /**
     Constructor that reads in a vector that is assumed to be sorted
   */
  FindInterval(vec inx) {
    x = conv_to<stdvec>::from(inx);
  }
  int operator()(double xi, int previous = 0) {
    int index = std::lower_bound(x.begin()+previous, x.end(), xi) - x.begin();
    return (x[index]==xi) ? index+1 : index;
  }
  Iter find(double xi, int previous = 0) {
    return std::find(x.begin()+previous, x.end(), xi);
  }
  bool member(double xi, int previous = 0) {
    return find(xi,previous) != x.end();
  }
  bool member(Iter iter) {
    return iter != x.end();
  }
  int index(Iter iter) {
    return iter - x.begin();
  }
  int index(double xi, int previous = 0) {
    return find(xi,previous) - x.begin();
  }
private:
  stdvec x;
};
  
/**
   Construct block diagonal matrix
   @param bag An ordered bag of matrices to be used to form the block diagonal matrix
 */
template<class Type>
Mat<Type> bdiag(field<Mat<Type> > bag) {
  int nr=0, nc=0;
  for (size_t i=0; i<bag.n_elem; i++) {
    nr += bag(i).n_rows;
    nc += bag(i).n_cols;
  }
  Mat<Type> out(nr,nc);
  out.zeros();
  nr=0; nc=0;
  for (size_t i=0; i<bag.n_elem; i++) {
    out(span(nr,nr+bag(i).n_rows-1), span(nc,nc+bag(i).n_cols-1)) = bag(i);
    nr += bag(i).n_rows;
    nc += bag(i).n_cols;
  }
  return out;
}
/**
   Construct block diagonal matrix
   @param m0 First matrix
   @param m1 Second matrix
 */
template<class Type>
Mat<Type> bdiag(Mat<Type> m0, Mat<Type> m1) {
  field<Mat<Type> > bag(2);
  bag(0) = m0;
  bag(1) = m1;
  return bdiag(bag);
}

std::function<mat (vec)> Fprob(int nStates, imat indices) {
  return [nStates,indices](vec X) -> mat {
    mat out(nStates,indices.n_rows);
    out.zeros();
    for (size_t i=0; i<indices.n_rows; i++) {
      int from = indices(i,0);
      int to = indices(i,1);
      out(to,i) = X(from);
      out(from,i) = -X(from);
    }
    return out;
  };
}

std::function<cube (vec)> Fjac(int nStates, std::function<mat(vec)> F) {
  return [nStates,F](vec X) -> cube {
    mat f = F(linspace(1.0,double(nStates),nStates));
    cube out(f.n_rows,f.n_rows,f.n_cols);
    for (size_t j=0; j<f.n_rows; j++) {
      vec ej(nStates);
      ej.zeros();
      ej(j) = 1.0;
      mat Fj = F(ej);
      for (size_t i=0; i<f.n_cols; i++) {
	for (size_t k=0; k<f.n_rows; k++)
	  out(k,j,i) = Fj(k,i);
      }
    }
    return out;
  };
}

std::function<mat (vec)> Fcombined(int nObs, int nStates, std::function<mat (vec)> F) {
  return [nObs,nStates,F](vec X) -> mat {
    field<mat> set(nObs);
    for (int i=0; i<nObs; i++) {
      set(i) = F(X(span(nStates*i,nStates*(i+1)-1)));
    }
    return bdiag<double>(set);
  };
}

std::function<mat(vec)> addFlos(std::function<mat (vec)> F) {
  return [F](vec X) -> mat {
    vec X1 = X(span(0,X.n_elem/2-1));
    mat F1 = F(X1);
    mat F2(X.n_elem/2,1);
    F2 = X1;
    return bdiag(F1,F2);
  };
}

mat makeW(int nStates, vec weights, int nOuter=1) {
  int nObs = weights.n_elem;
  mat W=eye(nStates,nStates)*weights(0);
  if (nObs>1)
    for (int i=1; i<nObs; i++)
      W = join_cols(W,eye(nStates,nStates)*weights(i));
  if (nOuter>1) {
    mat Wi = W;
    for (int i=1; i<nOuter; i++)
      W = bdiag(W,Wi);
  }
  return W;
}

/**
   Ryalen's plugin estimator using stochastic differential equations
   Discrete form
   @param n Number of observations in the initial analysis dataset
   @param hazMatrix Hazards for each transition (by row) for each time-point (by column)
   @param f A function that takes the state vector X and returns a matrix that is multiplied by the cum. hazard steps 
            at each time step to calculate the update in X
   @param gradf A function that takes the state vector X and returns a cube that is multiplied by the cum. hazard steps 
                at each time step to calculate the change in gradient for X
   @param X0 Initial values for the state vector
   @param V0 Initial values for the variance-covariance matrix (usually zero(nStates,nStates))
   @param W Weight matrix used to calculate summary measures
   @param vcov Boolean for whether to return the full variance-covariance matrix -- default=false, as this can be large
 */
PluginEstimateDiscrete
pluginEstimateDiscrete(int n, mat hazMatrix,
		       std::function<mat (vec)> f, std::function<cube (vec)> gradf,
		       vec X0, mat V0, mat W=zeros(1,1), bool vcov=false) {
  int numIncrements = hazMatrix.n_cols;
  mat X = zeros(X0.n_elem, numIncrements);
  cube covariance = vcov ? zeros(V0.n_rows, V0.n_cols, numIncrements) : zeros(1,1,1);
  mat variance = zeros(V0.n_rows, numIncrements);
  vec Xn = X0;
  mat Vn = V0;
  mat Y = zeros(1,1); // default if W is not specified
  mat varY = zeros(1,1); // default if W is not specified
  if (W.n_rows == X0.n_elem) {
    Y = zeros(W.n_cols,numIncrements);
    varY = zeros(W.n_cols,numIncrements);
  }
  cube gradf0 = gradf(X0); // for dimensions
  int Xrows = X0.n_elem,
    // hazRows = hazMatrix.n_rows,
    numHaz = gradf0.n_slices;
  X.col(0) = X0;
  if (vcov)
    covariance.slice(0) = V0;
  variance.col(0) = V0.diag();
  for (int i=1; i<numIncrements; i++) {
    mat fx = f(Xn);
    vec Xnew = Xn + fx * hazMatrix.col(i);
    mat Vtemp = zeros(Xrows, Xrows);
    cube fjac = gradf(Xn); // (a,b,c) for state-row a wrt state-var b in rate-col c
    for (int j=0; j<numHaz; j++) {
      Vtemp += (Vn*fjac.slice(j).t() + fjac.slice(j) * Vn) * hazMatrix(j,i);
    }
    vec dB = hazMatrix.col(i);
    mat dW = dB * dB.t();
    mat Vnew = Vn + Vtemp + n*fx*dW*fx.t();
    X.col(i) = Xnew;
    // if (i==1) {
    //   Rprintf("W.n_rows = %i\n", W.n_rows);
    //   Rprintf("W.n_cols = %i\n", W.n_cols);
    //   Rprintf("X0.n_elem = %i\n", X0.n_elem);
    // }
    if (W.n_rows == X0.n_elem) {
      Y.col(i) = W.t() * Xnew;
      varY.col(i) = mat(W.t() * Vnew * W).diag();
    }
    if (vcov)
      covariance.slice(i) = Vnew;
    variance.col(i) = Vnew.diag();
    Xn = Xnew;
    Vn = Vnew;
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
  }
  PluginEstimateDiscrete out = {X,variance/n,covariance/n,vcov,n,Y,varY};
  return out;
}

// TODO: use one function for both discrete and continuous case

/**
   Ryalen's plugin estimator using stochastic differential equations
   Mixed continuous/discrete form
   @param n Number of observations in the initial analysis dataset
   @param hazMatrix Hazards for each transition (by row) for each time-point (by column)
   @param f A function that takes the state vector X and returns a matrix that is multiplied by the cum. hazard steps 
            at each time step to calculate the update in X
   @param gradf A function that takes the state vector X and returns a cube that is multiplied by the cum. hazard steps 
                at each time step to calculate the change in gradient for X
   @param X0 Initial values for the state vector
   @param V0 Initial values for the variance-covariance matrix (usually zero(nStates,nStates))
   @param times Vector of the times, including the event times and a grid of evaluation times
   @param nOut Number of grid evaluation points to be output
   @param W Weight matrix used to calculate summary measures
   @param vcov Boolean for whether to return the full variance-covariance matrix -- default=false, as this can be large
   @param nLebesgue Number of full grid evaluation points
*/
PluginEstimateCts
pluginEstimateCts(int n, mat hazMatrix, std::function<mat(vec)> f, std::function<cube(vec)> gradf,
		vec X0, mat V0, vec times, int nOut=300, mat W = zeros(1,1), bool vcov=false, int nLebesgue=10001) {
  double // start=min(times),
    finish=max(times);
  // currently assumes start=0
  vec s = linspace(0,finish,nLebesgue);
  uvec sIndex = linspace<uvec>(0,nLebesgue-1,nOut);
  vec allTimes = unique(join_cols(s,times));
  vec subTimes = unique(join_cols(s(sIndex),times));
  FindInterval find_time(times); // assumes hazMatrix.n_cols == times.n_elem
  FindInterval find_subtime(subTimes);
  vec ds = diff(join_cols(vec({0.0}),allTimes));
  int numIncrements = subTimes.n_elem;
  mat X = zeros(X0.n_elem, numIncrements);
  cube covariance = vcov ? zeros(V0.n_rows, V0.n_cols, numIncrements) : zeros(1,1,1);
  mat variance = zeros(V0.n_rows, numIncrements);
  mat Y = zeros(1,1); // default if W is not specified
  mat varY = zeros(1,1); // default if W is not specified
  if (W.n_rows == X0.n_elem) {
    Y = zeros(W.n_cols,numIncrements);
    varY = zeros(W.n_cols,numIncrements);
  }
  vec Xn = X0;
  mat Vn = V0;
  cube gradf0 = gradf(X0); // for dimensions
  int Xrows = X0.n_elem,
    hazRows = hazMatrix.n_rows,
    numHaz = gradf0.n_slices;
  // bool Lebesgue = numHaz == hazRows+1; // ASSUMED TO BE TRUE
  int // k=0, // index in s
    l=1; // index in output
  X.col(0) = X0;
  if (vcov)
    covariance.slice(0) = V0;
  variance.col(0) = V0.diag();
  for (size_t i=1; i<allTimes.n_elem; i++) {
    double time = allTimes(i);
    FindInterval::Iter iter = find_time.find(time);
    bool event = find_time.member(iter);
    vec dsi = {ds(i)};
    vec hazVec;
    if (event) {
      vec haz = hazMatrix.col(find_time.index(iter));
      hazVec = join_cols(haz, dsi);
    } else {
      hazVec = join_cols(zeros(hazRows), dsi);
    }
    mat fx = f(Xn);
    vec Xnew = Xn + fx * hazVec;
    mat Vtemp = zeros(Xrows, Xrows);
    cube fjac = gradf(Xn); // (a,b,c) for state-row a wrt state-var b in rate-col c
    for (int j=0; j<numHaz; j++) {
      Vtemp += (Vn*fjac.slice(j).t() + fjac.slice(j) * Vn) * hazVec(j);
    }
    vec dB = hazVec;
    dB(dB.n_elem-1) = 0.0;
    mat dW = dB * dB.t();
    mat Vnew = Vn + Vtemp + n*fx*dW*fx.t();
    if (find_subtime.member(time)) {
      X.col(l) = Xnew;
      if (vcov)
	covariance.slice(l) = Vnew;
      variance.col(l) = Vnew.diag();
      if (W.n_rows == X0.n_elem) {
	Y.col(l) = W.t() * Xnew;
	varY.col(l) = mat(W.t() * Vnew * W).diag();
      }
      l++;
    }
    Xn = Xnew;
    Vn = Vnew;
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
  }
  PluginEstimateCts out = {X,variance/n,covariance/n,vcov,n,subTimes,Y,varY};
  return out;
}

// PluginEstimateDiscrete
// calc_P_by(int n, int nNewdata, mat hazMatrix, 
// 	  vec X0, imat tmat, vec weights, bool vcov=false) {
RcppExport SEXP plugin_P_by(SEXP _n, SEXP _nNewdata, SEXP _hazMatrix, 
			  SEXP _X0, SEXP _tmat, SEXP _weights, SEXP _vcov) {
  int n = Rcpp::as<int>(_n);
  size_t nNewdata = Rcpp::as<size_t>(_nNewdata);
  mat hazMatrix = Rcpp::as<mat>(_hazMatrix);
  vec X0 = Rcpp::as<vec>(_X0);
  imat tmat = Rcpp::as<imat>(_tmat);
  vec weights = Rcpp::as<vec>(_weights);
  bool vcov = Rcpp::as<bool>(_vcov);
  size_t nStates = tmat.max() - tmat.min() + 1; // assumes tmat is well formed (fragile)
  if (nStates == X0.n_elem) X0 = vec(repmat(mat(X0),nNewdata,1));
  std::function<mat(vec)> F = Fcombined(nNewdata,nStates,Fprob(nStates, tmat));
  std::function<cube(vec)> Fgrad = Fjac(nStates*nNewdata,F);
  mat V0 = zeros(nStates*nNewdata,nStates*nNewdata);
  mat W = (weights.n_elem == nNewdata) ? makeW(nStates, weights) : zeros(1,1);
  return Rcpp::wrap(pluginEstimateDiscrete(n, hazMatrix, F, Fgrad, X0, V0, W, vcov));
}

// PluginEstimateCts
// calc_P_L_by(int n, int nNewdata, mat hazMatrix, 
// 	    vec X0, imat tmat, vec time, vec weights, int nOut=300, bool vcov=false, int nLebesgue=10001) {
RcppExport SEXP plugin_P_L_by(SEXP _n, SEXP _nNewdata, SEXP _hazMatrix, 
			    SEXP _X0, SEXP _tmat, SEXP _time, SEXP _weights, SEXP _nOut,
			    SEXP _vcov, SEXP _nLebesgue) {
  int n = Rcpp::as<int>(_n);
  size_t nNewdata = Rcpp::as<size_t>(_nNewdata);
  mat hazMatrix = Rcpp::as<mat>(_hazMatrix);
  vec X0 = Rcpp::as<vec>(_X0);
  imat tmat = Rcpp::as<imat>(_tmat);
  vec time = Rcpp::as<vec>(_time);
  vec weights = Rcpp::as<vec>(_weights);
  int nOut = Rcpp::as<int>(_nOut);
  bool vcov = Rcpp::as<bool>(_vcov);
  int nLebesgue = Rcpp::as<int>(_nLebesgue);
  size_t nStates = tmat.max() - tmat.min() + 1; // assumes tmat is well formed (fragile)
  if (nStates == X0.n_elem)
    X0 = vec(repmat(mat(X0),nNewdata,1));
  vec X1 = join_cols(X0,zeros(X0.n_elem));
  std::function<mat(vec)> F = addFlos(Fcombined(nNewdata,nStates,Fprob(nStates, tmat)));
  std::function<cube(vec)> Fgrad = Fjac(nStates*nNewdata*2,F);
  mat V1 = zeros(nStates*nNewdata*2,nStates*nNewdata*2);
  mat W = (weights.n_elem == nNewdata) ? makeW(nStates, weights, 2) : zeros(1,1);
  return Rcpp::wrap(pluginEstimateCts(n, hazMatrix, F, Fgrad, X1, V1, time, nOut, W, vcov, nLebesgue));
}

// // [[Rcpp::depends(RcppArmadillo)]]
// //  [[Rcpp::export]]
// PluginEstimateDiscrete
// plugin_P(int n, mat hazMatrix, 
//       vec X0, imat tmat, bool vcov=false) {
//   int nStates = X0.n_elem;
//   std::function<mat(vec)> F = Fprob(nStates, tmat);
//   std::function<cube(vec)> Fgrad = Fjac(nStates,F);
//   mat V0 = zeros(nStates,nStates);
//   return pluginEstimateDiscrete(n, hazMatrix, F, Fgrad, X0, V0, zeros(1,1), vcov);
// }
// // [[Rcpp::depends(RcppArmadillo)]]
// //  [[Rcpp::export]]
// PluginEstimateCts
// plugin_P_L(int n, mat hazMatrix, 
// 	 vec X0, imat tmat, vec time, int nOut=300, bool vcov=false, int nLebesgue=10001) {
//   int m = X0.n_elem;
//   vec X1 = join_cols(X0,zeros(m));
//   std::function<mat(vec)> F = addFlos(Fprob(m, tmat));
//   std::function<cube(vec)> Fgrad = Fjac(2*m,F);
//   mat V1 = zeros(2*m,2*m);
//   return pluginEstimateCts(n, hazMatrix, F, Fgrad, X1, V1, time, nOut, zeros(1,1), vcov, nLebesgue);
// }
