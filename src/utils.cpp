#include <Rcpp.h>
#include <R_ext/Applic.h>

namespace rstpm2 {

  // Complete Q matrix from a QR decomposition
  Rcpp::NumericMatrix qr_q(const Rcpp::NumericMatrix& X, double tol) 
  {
    // Initialize member data and allocate heap memory
    int n=X.rows(), p=X.cols(), rank=0;
    Rcpp::NumericMatrix qr(X), y(n,n), q(n,n);
    int* pivot=(int*)R_alloc(p,sizeof(int)); 
    double* tau=(double*)R_alloc(p,sizeof(double)); 
    double* work=(double*)R_alloc(p*2,sizeof(double));
    for(int i=0;i<p;i++) 
      pivot[i]=i+1; 
    for(int i=0;i<n;i++) 
      for(int j=0;j<n;j++) 
	y(i,j) = i==j ? 1.0 : 0.0;
    // LINPACK QR factorization via householder transformations
    F77_CALL(dqrdc2)(&qr[0], &n, &n, &p, &tol, &rank, tau, pivot, work);
    // Compute orthogonal factor Q
    F77_CALL(dqrqy)(&qr[0], &n, &rank, tau, &y[0], &n, &q[0]);
    return q;
  }

}
