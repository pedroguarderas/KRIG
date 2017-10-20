
#include "RKHEstimate.h"


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
List RKHEstimate( const arma::mat& Z, const arma::mat& X, const arma::mat& Y,
                  const arma::mat& K, const arma::mat& k, const arma::mat& S,
                  const int type, const int cinv ) {

  int n = X.n_rows;
  int m = Y.n_rows;
  int p = Z.n_cols;
  
  arma::mat J( n, n );
  arma::mat W( m, p );

  if ( type == 2 ) {
    K = K + S;
  }
  
  if ( cinv == 0 ) {
    J = inv_sympd( K );
  } else if ( cinv == 1 ) {
    J = inv( K );
  } else if ( cinv == 2 ) {
    J = chol( K );
    J = inv( J );
    J = J.t() * J;
  }
  
  if ( type == 0 ) {
    W = k * J * Z;  
  } else if ( type == 1 ) {
    double v;
    arma::mat U( n, p );
    arma::mat u = arma::ones( n, 1 );
    arma::mat w = arma::ones( m, 1 );
    
    v = 1.0 / as_scalar( u.t() * J * u );
    U = v * u.t() * J * Z;
    W = w * U + k * J * ( Z - u * U );
  } else if ( type == 2 ) {
    double v;
    arma::mat U( n, p );
    arma::mat u = arma::ones( n, 1 );
    arma::mat w = arma::ones( m, 1 );
    
    v = 1.0 / as_scalar( u.t() * J * u );
    U = v * u.t() * J * Z;
    W = w * U + k * J * ( Z - u * U );
  }
  
  return List::create( Named( "W" ) = W,
                       Named( "K" ) = K,
                       Named( "k" ) = k,
                       Named( "J" ) = J );
}
