
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
List RKHEstimate( const arma::mat& Z, const arma::mat& K, const arma::mat& k,
                  const arma::mat& G, 
                  const arma::mat& g,
                  const std::string type, 
                  const std::string typeinv ) {

  int n = Z.n_rows;
  int m = k.n_cols;

  List KRIG;
  
  arma::mat J( n, n );
  arma::mat W( m, n );

  // Inverse computation
  if ( typeinv == "syminv" ) {
    J = inv_sympd( K );
    
  } else if ( typeinv == "inv" ) {
    J = inv( K );
    
  } else if ( typeinv == "cholinv" ) {
    J = chol( K );
    J = inv( J );
    J = J.t() * J;
    
  }
  
  // Kriging computation
  if ( type == "simple" ) { // Simple kriging
    W = k.t() * J * Z;  
    
    KRIG[ "Z" ] = W;
    KRIG[ "J" ] = J;
    
  } else if ( type == "ordinary" ) {  // Ordinary kriging
    double alpha;
    arma::mat u = arma::ones( n, 1 );
    arma::mat tau( n, m );

    alpha = 1.0 / as_scalar( u.t() * J * u );
    tau = arma::ones( n, n ) * J * k - arma::ones( n, m );
    W = ( k.t() - alpha * tau.t() ) * J * Z;
    
    KRIG[ "Z" ] = W;
    KRIG[ "J" ] = J;
    KRIG[ "alpha" ] = alpha;
    KRIG[ "tau" ] = tau;
    
  } else if ( type == "universal" ) { // Universal kriging
    
    int p = G.n_rows;
    arma::mat A( n, p );
    arma::mat tau( p, 1 );

    A = G.t() * inv_sympd( G * J * G.t() );
    tau = G.t() * K * k - g;
    W = ( k.t() - tau.t() * A.t() ) * J * Z;
    
    KRIG[ "Z" ] = W;
    KRIG[ "J" ] = J;
    KRIG[ "A" ] = A;
    KRIG[ "tau" ] = tau;

  }
  
  return KRIG;

}

