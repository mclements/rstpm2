#include <RcppArmadillo.h>
#include <errno.h>
#include <math.h>
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <algorithm>

/*
 *	Natural Splines
 *	---------------
 *	Here the end-conditions are determined by setting the second
 *	derivative of the spline at the end-points to equal to zero.
 *
 *	There are n-2 unknowns (y[i]'' at x[2], ..., x[n-1]) and n-2
 *	equations to determine them.  Either Choleski or Gaussian
 *	elimination could be used.
 */
static void
natural_spline(int n, double *x, double *y, double *b, double *c, double *d)
{
    int nm1, i;
    double t;
    x--; y--; b--; c--; d--;
    if(n < 2) {
	errno = EDOM;
	return;
    }
    if(n < 3) {
	t = (y[2] - y[1]);
	b[1] = t / (x[2]-x[1]);
	b[2] = b[1];
	c[1] = c[2] = d[1] = d[2] = 0.0;
	return;
    }
    nm1 = n - 1;
    /* Set up the tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */
    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1])/d[1];
    for( i=2 ; i<n ; i++) {
	d[i] = x[i+1] - x[i];
	b[i] = 2.0 * (d[i-1] + d[i]);
	c[i+1] = (y[i+1] - y[i])/d[i];
	c[i] = c[i+1] - c[i];
    }
    /* Gaussian elimination */
    for(i=3 ; i<n ; i++) {
	t = d[i-1]/b[i-1];
	b[i] = b[i] - t*d[i-1];
	c[i] = c[i] - t*c[i-1];
    }
    /* Backward substitution */
    c[nm1] = c[nm1]/b[nm1];
    for(i=n-2 ; i>1 ; i--)
	c[i] = (c[i]-d[i]*c[i+1])/b[i];
    /* End conditions */
    c[1] = c[n] = 0.0;
    /* Get cubic coefficients */
    b[1] = (y[2] - y[1])/d[1] - d[i] * c[2];
    c[1] = 0.0;
    d[1] = c[2]/d[1];
    b[n] = (y[n] - y[nm1])/d[nm1] + d[nm1] * c[nm1];
    for(i=2 ; i<n ; i++) {
	b[i] = (y[i+1]-y[i])/d[i] - d[i]*(c[i+1]+2.0*c[i]);
	d[i] = (c[i+1]-c[i])/d[i];
	c[i] = 3.0*c[i];
    }
    c[n] = 0.0;
    d[n] = 0.0;
    return;
}
/*
 *	Splines a la Forsythe Malcolm and Moler
 *	---------------------------------------
 *	In this case the end-conditions are determined by fitting
 *	cubic polynomials to the first and last 4 points and matching
 *	the third derivitives of the spline at the end-points to the
 *	third derivatives of these cubics at the end-points.
 */
static void
fmm_spline(int n, double *x, double *y, double *b, double *c, double *d)
{
    int nm1, i;
    double t;
    /* Adjustment for 1-based arrays */
    x--; y--; b--; c--; d--;
    if(n < 2) {
	errno = EDOM;
	return;
    }
    if(n < 3) {
	t = (y[2] - y[1]);
	b[1] = t / (x[2]-x[1]);
	b[2] = b[1];
	c[1] = c[2] = d[1] = d[2] = 0.0;
	return;
    }
    nm1 = n - 1;
    /* Set up tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */
    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1])/d[1];/* = +/- Inf	for x[1]=x[2] -- problem? */
    for(i=2 ; i<n ; i++) {
	d[i] = x[i+1] - x[i];
	b[i] = 2.0 * (d[i-1] + d[i]);
	c[i+1] = (y[i+1] - y[i])/d[i];
	c[i] = c[i+1] - c[i];
    }
    /* End conditions. */
    /* Third derivatives at x[0] and x[n-1] obtained */
    /* from divided differences */
    b[1] = -d[1];
    b[n] = -d[nm1];
    c[1] = c[n] = 0.0;
    if(n > 3) {
	c[1] = c[3]/(x[4]-x[2]) - c[2]/(x[3]-x[1]);
	c[n] = c[nm1]/(x[n] - x[n-2]) - c[n-2]/(x[nm1]-x[n-3]);
	c[1] = c[1]*d[1]*d[1]/(x[4]-x[1]);
	c[n] = -c[n]*d[nm1]*d[nm1]/(x[n]-x[n-3]);
    }
    /* Gaussian elimination */
    for(i=2 ; i<=n ; i++) {
	t = d[i-1]/b[i-1];
	b[i] = b[i] - t*d[i-1];
	c[i] = c[i] - t*c[i-1];
    }
    /* Backward substitution */
    c[n] = c[n]/b[n];
    for(i=nm1 ; i>=1 ; i--)
	c[i] = (c[i]-d[i]*c[i+1])/b[i];
    /* c[i] is now the sigma[i-1] of the text */
    /* Compute polynomial coefficients */
    b[n] = (y[n] - y[n-1])/d[n-1] + d[n-1]*(c[n-1]+ 2.0*c[n]);
    for(i=1 ; i<=nm1 ; i++) {
	b[i] = (y[i+1]-y[i])/d[i] - d[i]*(c[i+1]+2.0*c[i]);
	d[i] = (c[i+1]-c[i])/d[i];
	c[i] = 3.0*c[i];
    }
    c[n] = 3.0*c[n];
    d[n] = d[nm1];
    return;
}
/*
 *	Periodic Spline
 *	---------------
 *	The end conditions here match spline (and its derivatives)
 *	at x[1] and x[n].
 *
 *	Note: There is an explicit check that the user has supplied
 *	data with y[1] equal to y[n].
 */
static void
periodic_spline(int n, double *x, double *y, double *b, double *c, double *d)
{
    double s;
    int i, nm1;
    // double *e = (double *) R_alloc(n, sizeof(double));
    arma::vec ev(n);
    double* e = ev.memptr();
    /* Adjustment for 1-based arrays */
    x--; y--; b--; c--; d--; e--;
    if(n < 2 || y[1] != y[n]) {
	errno = EDOM;
	return;
    }
    if(n == 2) {
	b[1] = b[2] = c[1] = c[2] = d[1] = d[2] = 0.0;
	return;
    } else if (n == 3) {
	b[1] = b[2] = b[3] = -(y[1] - y[2])*(x[1] - 2*x[2] + x[3])/(x[3]-x[2])/(x[2]-x[1]);
	c[1] = -3*(y[1]-y[2])/(x[3]-x[2])/(x[2]-x[1]);
	c[2] = -c[1];
	c[3] = c[1];
	d[1] = -2*c[1]/3/(x[2]-x[1]);
	d[2] = -d[1]*(x[2]-x[1])/(x[3]-x[2]);
	d[3] = d[1];
	return;
    }
    /* else --------- n >= 4 --------- */
    nm1 = n-1;
    /* Set up the matrix system */
    /* A = diagonal	 B = off-diagonal  C = rhs */
#define A	b
#define B	d
#define C	c
    B[1]  = x[2] - x[1];
    B[nm1]= x[n] - x[nm1];
    A[1] = 2.0 * (B[1] + B[nm1]);
    C[1] = (y[2] - y[1])/B[1] - (y[n] - y[nm1])/B[nm1];
    for(i = 2; i < n; i++) {
	B[i] = x[i+1] - x[i];
	A[i] = 2.0 * (B[i] + B[i-1]);
	C[i] = (y[i+1] - y[i])/B[i] - (y[i] - y[i-1])/B[i-1];
    }
    /* Choleski decomposition */
#define L	b
#define M	d
#define E	e
    L[1] = sqrt(A[1]);
    E[1] = (x[n] - x[nm1])/L[1];
    s = 0.0;
    for(i = 1; i <= nm1 - 2; i++) {
	M[i] = B[i]/L[i];
	if(i != 1) E[i] = -E[i-1] * M[i-1] / L[i];
	L[i+1] = sqrt(A[i+1]-M[i]*M[i]);
	s = s + E[i] * E[i];
    }
    M[nm1-1] = (B[nm1-1] - E[nm1-2] * M[nm1-2])/L[nm1-1];
    L[nm1] = sqrt(A[nm1] - M[nm1-1]*M[nm1-1] - s);
    /* Forward Elimination */
#define Y	c
#define D	c
    Y[1] = D[1]/L[1];
    s = 0.0;
    for(i=2 ; i<=nm1-1 ; i++) {
	Y[i] = (D[i] - M[i-1]*Y[i-1])/L[i];
	s = s + E[i-1] * Y[i-1];
    }
    Y[nm1] = (D[nm1] - M[nm1-1] * Y[nm1-1] - s) / L[nm1];
#define X	c
    X[nm1] = Y[nm1]/L[nm1];
    X[nm1-1] = (Y[nm1-1] - M[nm1-1] * X[nm1])/L[nm1-1];
    for(i=nm1-2 ; i>=1 ; i--)
	X[i] = (Y[i] - M[i] * X[i+1] - E[i] * X[nm1])/L[i];
    /* Wrap around */
    X[n] = X[1];
    /* Compute polynomial coefficients */
    for(i=1 ; i<=nm1 ; i++) {
	s = x[i+1] - x[i];
	b[i] = (y[i+1]-y[i])/s - s*(c[i+1]+2.0*c[i]);
	d[i] = (c[i+1]-c[i])/s;
	c[i] = 3.0*c[i];
    }
    b[n] = b[1];
    c[n] = c[1];
    d[n] = d[1];
    return;
}
#undef A
#undef B
#undef C
#undef L
#undef M
#undef E
#undef Y
#undef D
#undef X
/* These were the public interfaces */
static void
spline_coef(int method, int n, double *x, double *y,
	    double *b, double *c, double *d)
{
    switch(method) {
    case 1:
    	periodic_spline(n, x, y, b, c, d);	break;
    case 2:
	natural_spline(n, x, y, b, c, d);	break;
    case 3:
	fmm_spline(n, x, y, b, c, d);	break;
    }
}
static void
spline_eval(int method, int nu, double *u, double *v,
	    int n, double *x, double *y, double *b, double *c, double *d)
{
/* Evaluate  v[l] := spline(u[l], ...),	    l = 1,..,nu, i.e. 0:(nu-1)
 * Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(*n-1)
 */
    const int n_1 = n - 1;
    int i, j, k, l;
    double ul, dx, tmp;
    if(method == 1 && n > 1) { /* periodic */
    	dx = x[n_1] - x[0];
    	for(l = 0; l < nu; l++) {
    	    v[l] = fmod(u[l]-x[0], dx);
    	    if(v[l] < 0.0) v[l] += dx;
    	    v[l] += x[0];
    	}
    } else
      for(l = 0; l < nu; l++) v[l] = u[l];
    for(l = 0, i = 0; l < nu; l++) {
	ul = v[l];
	if(ul < x[i] || (i < n_1 && x[i+1] < ul)) {
	    /* reset i  such that  x[i] <= ul <= x[i+1] : */
	    i = 0;
	    j = n;
	    do {
		k = (i+j)/2;
		if(ul < x[k]) j = k;
		else i = k;
	    } while(j > i+1);
	}
	dx = ul - x[i];
	/* for natural splines extrapolate linearly left */
	tmp = (method == 2 && ul < x[0]) ? 0.0 : d[i];
	v[l] = y[i] + dx*(b[i] + dx*(c[i] + dx*tmp));
    }
}
class SplineCoef {
public:
  SplineCoef(arma::vec x, arma::vec y, int method = 2) : x(x), y(y), method(method) {
    int n = x.n_elem;
    // if(y.n_elem != x.n_elem) REprintf("inputs of different lengths");
    b.resize(n); c.resize(n); d.resize(n);
    b.zeros(); c.zeros(); d.zeros();
    if (method==1 && y(0) != y(n-1))
      y(n-1)=y(0);
    spline_coef(method, x.n_elem, x.memptr(), y.memptr(), b.memptr(), c.memptr(), d.memptr());
  }
  arma::vec eval(arma::vec xout) {
    int nu = xout.n_elem;
    arma::vec yout(nu);
    spline_eval(method, nu, xout.memptr(), yout.memptr(),
		x.n_elem, x.memptr(), y.memptr(), b.memptr(), c.memptr(), d.memptr());
    return yout;
  }
  double eval(double xout) {
    double yout;
    spline_eval(method, 1, &xout, &yout,
		x.n_elem, x.memptr(), y.memptr(), b.memptr(), c.memptr(), d.memptr());
    return yout;
  }
  static arma::vec eval(double xout, std::vector<SplineCoef>& z);
  arma::vec x, y, b, c, d;
  int method;
};
arma::vec SplineCoef::eval(double xout, std::vector<SplineCoef>& z)
  {
    arma::vec yout(z.size());
    for (size_t i=0; i<z.size(); i++)
      yout(i) = z[i].eval(xout);
    return yout;
  }
namespace boost { namespace numeric { namespace odeint {
      template <>
      struct is_resizeable<arma::vec>
      {
	typedef boost::true_type type;
	const static bool value = type::value;
      };
      template <>
      struct same_size_impl<arma::vec, arma::vec>
      {
	static bool same_size(const arma::vec& x, const arma::vec& y)
	{
	  return x.size() == y.size();
	}
      };
      template<>
      struct resize_impl<arma::vec, arma::vec>
      {
	static void resize(arma::vec& v1, const arma::vec& v2)
	{
	  v1.resize(v2.size());
	}
      };
    } } } // namespace boost::numeric::odeint
typedef arma::vec state_type;
typedef std::vector<state_type> vector_state_type;
typedef std::vector<double> vector_times;
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }
    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
struct push_back_state
{
    std::vector< state_type >& m_states;
    push_back_state( std::vector< state_type > &states)
    : m_states( states ) { }
    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
    }
};
struct Flow {
  size_t from, to;
  SplineCoef s;
  std::vector<SplineCoef> gradients;
  bool use_log;
};
// old version of libarmadillo-dev: mat.as_col() not available
arma::vec flatten(arma::mat m) {
  arma::vec v(m.n_cols*m.n_rows);
  for (size_t j=0; j<m.n_cols; j++)
    for (size_t i=0; i<m.n_rows; i++)
      v(j*m.n_rows+i) = m(i,j);
  return v;
}
struct StateComponents {
  arma::vec P, L;
  arma::mat gradP, gradL;
};
struct StateComponentsCombined {
  arma::vec times;
  arma::mat P, L;
  arma::cube gradP, gradL;
};
// forward declaration
class Markov {
public:
  typedef std::vector<Flow> Flows;
  Flows& flows;
  double minTm;
  size_t nStates, nGradients;
  StateComponents stateComponents(const state_type& x) {
    arma::vec P = x(arma::span(0,nStates-1));
    arma::vec L = x(arma::span(nStates,2*nStates-1));
    arma::mat gradP(arma::vec(x(arma::span(2*nStates,2*nStates+nStates*nGradients-1))).memptr(),
		    nStates, nGradients);
    arma::mat gradL(arma::vec(x(arma::span(2*nStates+nStates*nGradients,
					   2*nStates+2*nStates*nGradients-1))).memptr(),
		    nStates, nGradients);
    return {P, L, gradP, gradL};
  }
  Markov(Flows& flows, double minTm = 1.0e-8) : flows(flows), minTm(minTm) {
    nGradients = nStates = 0;
    for(size_t i=0; i<flows.size(); i++) {
      nGradients += flows[i].gradients.size();
      nStates = std::max(std::max(nStates,flows[i].from),flows[i].to);
    }
    nStates++;
  }
  void operator() ( const state_type &x , state_type &dxdt , const double t )
  {
    arma::vec rates(flows.size());
    arma::vec dPdt(nStates,arma::fill::zeros);
    arma::vec dLdt;
    arma::mat dGradPdt(nStates,nGradients,arma::fill::zeros);
    arma::mat dGradLdt;
    // destructure the state vector
    StateComponents c = stateComponents(x);
    // P (transition probabilities)
    for (size_t i=0; i<flows.size(); i++) {
      rates(i) = flows[i].s.eval(flows[i].use_log ? log(t<minTm ? minTm : t) : t);
      // if (flows[i].use_log) rates(i) = exp(rates(i));
      double delta = c.P(flows[i].from)*rates(i);
      dPdt(flows[i].to) += delta;
      dPdt(flows[i].from) -= delta;
    }
    // L (length of stay)
    dLdt = c.P;
    // gradients for P
    for (size_t i=0, start=0; i<flows.size(); start+=flows[i].gradients.size(), i++) {
      for (size_t j=0; j<flows[i].gradients.size(); j++) {
	double delta = c.gradP(flows[i].from,start+j)*rates(i)+
	  c.P(flows[i].from)*flows[i].gradients[j].eval(t);
	dGradPdt(flows[i].to,start+j) += delta;
	dGradPdt(flows[i].from,start+j) -= delta;
      }
    }
    // gradients for L
    dGradLdt = c.gradP;
    dxdt = arma::join_cols(dPdt,
			   arma::join_cols(dLdt,
					   arma::join_cols(flatten(dGradPdt),
							   flatten(dGradLdt))));
  }
};
StateComponentsCombined run(Markov &model, arma::vec p0, arma::vec times) {
  using namespace boost::numeric::odeint;
  vector_state_type states;
  // times, outtimes and vtimes all have the same data
  vector_times outtimes; // not strictly needed - we could use push_back_state()
  vector_times vtimes = arma::conv_to< vector_times >::from(times);
  state_type x = join_cols(p0, arma::zeros(model.nStates+2*model.nStates*model.nGradients));
  BOOST_STATIC_ASSERT( is_resizeable<state_type>::value == true );
  integrate_times(make_dense_output( 1.0e-10 , 1.0e-10 , runge_kutta_dopri5< state_type >() ),
		  model, x,
		  vtimes.begin(),
		  vtimes.end(),
		  1.0,
		  push_back_state_and_time(states, outtimes));
  size_t nx = states[0].size(), nTimes = times.size();
  // combine the results
  arma::mat combined(nTimes,nx);
  for (size_t i=0; i<nTimes; i++)
    combined(arma::span(i,i),arma::span(0,nx-1)) = states[i].t();
  // destructure the combined matrix
  arma::mat P = combined(arma::span::all,arma::span(0,model.nStates-1));
  arma::mat L = combined(arma::span::all,arma::span(model.nStates,2*model.nStates-1));
  arma::cube gradP(arma::mat(combined(arma::span::all,
				      arma::span(2*model.nStates,2*model.nStates+model.nStates*model.nGradients-1))).memptr(),
		   nTimes, model.nStates, model.nGradients);
  arma::cube gradL(arma::mat(combined(arma::span::all,
			    arma::span(2*model.nStates+model.nStates*model.nGradients,
				       2*model.nStates+2*model.nStates*model.nGradients-1))).memptr(),
		   nTimes, model.nStates, model.nGradients);
  return {times, P, L, gradP, gradL};
}
// Rcpp::List runMarkovODE(arma::vec y0, arma::vec times, arma::vec tlam, arma::mat lam, Rcpp::List gradients, arma::ivec from, arma::ivec to, Rcpp::LogicalVector use_logs, double minTm = 1.0e-8) {
RcppExport SEXP runMarkovODE(SEXP _y0, SEXP _times, SEXP _tlam, SEXP _lam, SEXP _gradients, SEXP _from, SEXP _to, SEXP _use_logs, SEXP _minTm) {
  using namespace Rcpp;
  arma::vec y0 = as<arma::vec>(_y0);
  arma::vec times = as<arma::vec>(_times);
  arma::vec tlam = as<arma::vec>(_tlam);
  arma::mat lam = as<arma::mat>(_lam);
  List gradients = as<List>(_gradients);
  arma::ivec from = as<arma::ivec>(_from);
  arma::ivec to = as<arma::ivec>(_to);
  LogicalVector use_logs = as<LogicalVector>(_use_logs);
  double minTm = as<double>(_minTm);
  Markov::Flows flows;
  for (size_t i=0; i<from.size(); i++) {
    SplineCoef coef = SplineCoef(use_logs[i] ? log(tlam) : tlam, lam.col(i));
    std::vector<SplineCoef> vgradients;
    arma::mat g = Rcpp::as<arma::mat>(gradients[i]);
    for (int j=0; j<gradients.size(); j++)
      vgradients.push_back(SplineCoef(tlam, g.col(j)));
    flows.push_back({size_t(from(i)), size_t(to(i)), coef, vgradients,
                     bool(use_logs[i])});
  }
  Markov markov(flows, minTm);
  auto report = run(markov, y0, times);
  return wrap(List::create(Named("times")=report.times,
			   Named("P")=report.P,
			   Named("L")=report.L,
			   Named("gradP")=report.gradP,
			   Named("gradL")=report.gradL));
}

RcppExport SEXP multistate_ddt(SEXP P_, SEXP Pu_, SEXP Q_, SEXP Qu_) {
  using namespace arma;
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


class ExpM {
public:
  arma::mat Qmat;
  ExpM(arma::mat _Qmat) : Qmat(_Qmat) { }
  void operator() ( const state_type &x , state_type &dxdt , const double t )
  {
    dxdt = (x.t() * Qmat).t();
  }
};
RcppExport SEXP runExpM(SEXP _y0, SEXP _times, SEXP _Qmat) {
  using namespace Rcpp;
  arma::vec y0 = as<arma::vec>(_y0);
  arma::vec times = as<arma::vec>(_times);
  arma::mat Qmat = as<arma::mat>(_Qmat);
  ExpM model(Qmat);
  using namespace boost::numeric::odeint;
  vector_state_type states;
  // times, outtimes and vtimes all have the same data
  vector_times outtimes; // not strictly needed - we could use push_back_state()
  vector_times vtimes = arma::conv_to< vector_times >::from(times);
  BOOST_STATIC_ASSERT( is_resizeable<state_type>::value == true );
  integrate_times(make_dense_output( 1.0e-10 , 1.0e-10 , runge_kutta_dopri5< state_type >() ),
		  model, y0,
		  vtimes.begin(),
		  vtimes.end(),
		  1.0,
		  push_back_state_and_time(states, outtimes));
  size_t nx = states[0].size(), nTimes = times.size();
  // combine the results
  arma::mat combined(nTimes,nx);
  for (size_t i=0; i<nTimes; i++)
    combined(arma::span(i,i),arma::span(0,nx-1)) = states[i].t();
  return wrap(combined);
}
