
#include "krig_cokriging.h"

//--------------------------------------------------------------------------------------------------
List coKrig( const arma::mat& Z, 
             const arma::mat& K, 
             const arma::cube& k,
             const arma::mat& G, 
             const arma::mat& g,
             const std::string type = "ordinary", 
             const std::string cinv = "syminv" ) {
  
  int nq = K.n_rows;
  int n = Z.n_rows;
  int m = k.n_cols;
  int q = k.n_slices;
  
  List KRIG;
  
  arma::mat J( nq, nq );
  arma::mat W( m, n );
  arma::mat L( m, n );
  
  // Inverse computation of the covariance matrix
  if ( cinv == "syminv" ) {
    J = inv_sympd( K );
    
  } else if ( cinv == "inv" ) {
    J = inv( K );
    
  } else if ( cinv == "cholinv" ) {
    J = chol( K );
    J = inv( J );
    J = J.t() * J;
    
  } else if ( cinv == "ginv" ) {
    J = K;
  }
  
  // Kriging computation
  if ( type == "simple" ) { // Simple kriging
    
    L = J *  k;
    W = L.t() * Z;  
    
    KRIG[ "Z" ] = W;
    KRIG[ "L" ] = L;
    KRIG[ "J" ] = J;
    
    
  } else if ( type == "ordinary" ) {  // Ordinary kriging
    double alpha;
    arma::mat u = arma::ones( n, 1 );
    arma::mat tau( n, m );
    
    alpha = 1.0 / as_scalar( u.t() * J * u );
    tau = arma::ones( n, m ) - arma::ones( n, n ) * J * k;
    
    L = J * ( k + alpha * tau ) ;
    W = L.t() * Z;
    
    KRIG[ "Z" ] = W;
    KRIG[ "L" ] = L;
    KRIG[ "J" ] = J;
    KRIG[ "alpha" ] = alpha;
    KRIG[ "tau" ] = tau;
    
  } else if ( type == "universal" ) { // Universal kriging
    
    int p = G.n_rows;
    arma::mat A( n, p );
    arma::mat tau( p, 1 );
    
    A = G.t() * inv_sympd( G * J * G.t() );
    tau = g - G * J * k;
    
    L = J * ( k + A * tau ) ;
    W = L.t() * Z;
    
    KRIG[ "Z" ] = W;
    KRIG[ "L" ] = L;
    KRIG[ "J" ] = J;
    KRIG[ "A" ] = A;
    KRIG[ "tau" ] = tau;
    
  }
  
  return KRIG;
  
}
