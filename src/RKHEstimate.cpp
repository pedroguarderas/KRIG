
#include "RKHEstimate.h"

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
arma::mat RKHCov( const arma::mat& X, const arma::mat& Y, Function Kern, const bool symmetric ) {
  int i, j;
  int m = X.n_rows;
  int n = Y.n_rows;
  arma::mat K( m, n );
  arma::rowvec x, y;
  
  // Filling Gaussian process covariance matrix
  if ( symmetric ) {
    for ( i = 0; i < m; i++ ) { 
      for ( j = i; j < n; j++ ) {
        x = X.row( i );
        y = Y.row( j );
        K( i, j ) = as<double>( Kern( x, y ) );
        
        if ( j > i ) {
          K( j, i ) = K( i, j );
        }
        
      }
    }
    
  } else { 
    for ( i = 0; i < m; i++ ) { 
      for ( j = 0; j < n; j++ ) {
        x = X.row( i );
        y = Y.row( j );
        K( i, j ) = as<double>( Kern( x, y ) );
      }
    }
  }
  
  return K;

}

//--------------------------------------------------------------------------------------------------
List RKHEstimate( const arma::mat& Z, const arma::mat& X, const arma::mat& Y, Function Kern ) {

  int n = X.n_rows;
  int m = Y.n_rows;
  int p = Z.n_cols;
  
  arma::mat k( m, n );
  arma::mat K( n, n );
  arma::mat J( n, n );
  arma::mat W( m, p );
  arma::mat U( n, p );
  arma::mat O = arma::ones( n, 1 );
  
  K = RKHCov( X, X, Kern, true );
  k = RKHCov( Y, X, Kern );
  
  J = inv_sympd( K );
  // W = k * J * Z;
  
  U = ( 1.0 / ( O.t() * J * O ) ) * ( O.t() * J * Z );
  W = U + k * J * ( Z - U );

  return List::create( Named( "W" ) = W,
                       Named( "K" ) = K,
                       Named( "k" ) = k,
                       Named( "U" ) = U,
                       Named( "J" ) = J );
}
