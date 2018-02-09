
#include "krig_cokriging.h"

//--------------------------------------------------------------------------------------------------
List coKrig( const arma::mat& Z, 
             const arma::mat& K, 
             const arma::cube& k,
             const arma::mat& G, 
             const arma::mat& g,
             const std::string type = "ordinary", 
             const std::string cinv = "syminv" ) {
 
   // The dimensions considerations
   // dim( K ) = nq x nq
   // dim( k_r ) = nq x m x q
   // dim( l_r ) = dim( K * k_r )  = nq x m
   // dim( Z ) = n x q => dim( vect( Z ) ) = nq x 1
   // dim( X ) = m x d
   // dim( W ) = m x q
  int nq = K.n_rows;
  int n = Z.n_rows;
  int m = k.n_cols;
  int q = k.n_slices;
  
  int k;
  
  List KRIG;
  
  arma::mat J( nq, nq );
  arma::mat W( m, q );
  arma::cube L( nq, m, q );
  
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
    
    for ( k = 0; k < q; k++ ) {
      L.slice( k ) = J * k.slice( k );
    }
    
  } else if ( type == "ordinary" ) {  // Ordinary kriging
    double alpha;
    arma::mat u = arma::ones( n, 1 );
    arma::mat tau( n, m );
    
    alpha = 1.0 / as_scalar( u.t() * J * u );
    tau = arma::ones( n, m ) - arma::ones( n, n ) * J * k;
    
    L = J * ( k + alpha * tau ) ;
    for ( k = 0; k < q; k++ ) {
      L.slice( k ) = J * k.slice( k );
    }+
    
  } else if ( type == "universal" ) { // Universal kriging
    
    int p = G.n_rows;
    arma::mat A( n, p );
    arma::mat tau( p, 1 );
    
    A = G.t() * inv_sympd( G * J * G.t() );
    tau = g - G * J * k;
    for ( k = 0; k < q; k++ ) {
      L.slice( k ) = J * k.slice( k );
    }
    
  }
  
  W = L.t() * Z;
  
  return KRIG;
  
}
