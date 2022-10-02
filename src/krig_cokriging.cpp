
#include "krig_cokriging.h"

//--------------------------------------------------------------------------------------------------
List coKrig( const Eigen::MatrixXd& Z, 
             const Eigen::MatrixXd& K, 
             const Eigen::ArrayXd& k,
             const Eigen::MatrixXd& G, 
             const Eigen::ArrayXd& g,
             const std::string type, 
             const std::string cinv ) {
 
   // The dimensions considerations
   // dim( K ) = nq x nq
   // dim( k_r ) = nq x m x q
   // dim( l_r ) = dim( K * k_r )  = nq x m
   // dim( Z ) = n x q => dim( vect( Z ) ) = nq x 1
   // dim( X ) = m x d
   // dim( W ) = m x q
  int nq = K.rows();
  int n = Z.rows();
  int m = k.cols();
  int q = k.n_slices;
  
  int i, j, r, s;
  
  List KRIG;
  
  Eigen::MatrixXd J( nq, nq );
  Eigen::MatrixXd W( m, q );
  Eigen::ArrayXd L( nq, m, q );
  
  // Inverse computation of the covariance matrix
  if ( cinv == "syminv" ) {
    J = K.inverse();
    
  } else if ( cinv == "inv" ) {
    J = K.inverse();
    
  } else if ( cinv == "cholinv" ) {
    J = K.ldlt();

  } else if ( cinv == "ginv" ) {
    J = K;
  }
  
  // Kriging computation
  if ( type == "simple" ) { // Simple kriging
    
    for ( r = 0; r < q; r++ ) {
      L.slice( r ) = J * k.slice( r );
    }
    
  } else if ( type == "ordinary" ) {  // Ordinary kriging
    double alpha;
    Eigen::MatrixXd u = Eigen::MatrixXd::Ones( nq, 1 );
    Eigen::MatrixXd tau( nq, m );
    
    alpha = 1.0 / u.transpose() * J * u;
    
    for ( r = 0; r < q; r++ ) {
      tau = Eigen::MatrixXd::Ones( nq, m ) - Eigen::MatrixXd::Ones( nq, nq ) * J * k.slice( r );
      L.slice( r ) = J * ( k.slice( r ) + alpha * tau );
    }
    
  } else if ( type == "universal" ) { // Universal kriging
    
    int p = G.rows();
    Eigen::MatrixXd A( nq, p );
    Eigen::MatrixXd tau( p, m );
    
    A = G.transpose() * inv_sympd( G * J * G.transpose() );
    for ( r = 0; r < q; r++ ) {
      tau = g.slice( r ) - G * J * k.slice( r );
      L.slice( r ) = J * ( k.slice( r ) + A * tau );
    }
    
  }
  
  // for ( i = 0; i < n; i++ ) {
  //   for ( r = 0; r < q; r++ ) {
  //     for( j = 0; j < q; j++ ) {
  //       l( r, j, i ) = L( i + j * ( n - 1 ), s, r );
  //     }
  //   }
  // }
    
  for ( r = 0; r < q; r++ ) {
    W.col( r ) = L.slice( r ).transpose() * vectorise( Z, 1 );
  }
  
  KRIG[ "Z" ] = W;
  KRIG[ "L" ] = L;
  KRIG[ "J" ] = J;
  
  return KRIG;
  
}
