#ifndef __RKHIntegral__
#define __RKHIntegral__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/*__________________________________________________________________________________________________
  Numerical integration funtion
 */
// [[Rcpp::export]]
List integralRKH( List Kernels ) {

  List kernel;
  List I( Kernels.size() );
  NumericVector X;
  NumericVector A( Kernels.size(), 0.0 );

  int i, j, k, N;
  double x, y, a, b, h, n;
  
  for ( k = 0; k < Kernels.size(); k++ ) {
    kernel = Kernels[k];
    Function f = kernel["kernel"]; 
    X = kernel["grid"];
    NumericVector IK( X.size(), 0.0 );
    b = kernel["max"];
    a = kernel["min"];
    n = kernel["size"];
    N = (int)n;
    h = ( b - a ) / ( n - 1.0 ); 

    // #pragma omp parallel for shared( a, b, n, h, N, X, IK, f ) private( x, i, j ) schedule( static ) collapse(2)
    for ( j = 0; j < X.size(); j++ ) {
      for ( i = 0; i < N; i++ ) {
        x = a + ( ( double )i ) * h;
        IK[j] = IK[j] + as<double>( f( x, X[j] ) ) * h;
      }
    }
    
    // #pragma omp parallel for shared( a, k, h, N, A ) private( x, y, i, j ) schedule( static ) collapse(2)
    for ( i = 0; i < N; i++ ) {
      for ( j = 0; j < N; j++ ) {
        x = a + ( ( double )i ) * h;
        y = a + ( ( double )j ) * h;
        A[k] = A[k] + as<double>( f( x, y ) ) * h * h;
      }
    }
    
    I[k] = List::create( Named( "integral" ) = IK, Named( "alpha" ) = A[k] );
  }
  
  return I;
}

/*__________________________________________________________________________________________________
  Kernels evaluation
 */
// [[Rcpp::export]]
List evalKernRKH( List Kernels, List Integral ) {
  
  int I, N, i, k, l;
  double alpha;
  List kernel, integral;
  NumericVector X, Int;

  I = Kernels.size();
  kernel = Kernels[0];
  X = kernel["grid"];
  N = X.size();
  
  List Gamma( I );
  NumericMatrix KANOVA( N, N );

  fill( KANOVA.begin(), KANOVA.end(), 1.0 );
  
  for ( i = 0; i < I; i++ ) {
    
    kernel = Kernels[i];
    integral = Integral[i];
    Function f = kernel["kernel"];
    X = kernel["grid"];
    Int = integral["integral"];
    alpha = integral["alpha"];
    NumericMatrix G( N, N );

    for ( k = 0; k < N; k++ ) {
      for ( l = k; l < N; l++ ) {

        G(k,l) = as<double>( f( X[k], X[l] ) ) - Int[k] * Int[l] / alpha;
        KANOVA(k,l) = KANOVA(k,l) * ( 1.0 + G(k,l) );

        if ( l > k ) {
          G(l,k) = G(k,l);
          KANOVA(l,k) = KANOVA(k,l);
        }
      }
    }
    Gamma[i] = G;
  }
  return List::create( Named( "Gamma" ) = Gamma, Named( "KANOVA" ) = KANOVA );
}

