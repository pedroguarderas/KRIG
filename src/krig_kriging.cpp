
#include "krig_kriging.h"


//--------------------------------------------------------------------------------------------------
Eigen::MatrixXd Kov( const Eigen::MatrixXd& X, 
                     const Eigen::MatrixXd& Y, 
                     Function Kern, 
                     const bool symmetric ) {
  int i, j;
  int m = X.rows();
  int n = Y.rows();
  Eigen::MatrixXd K( n, m );
  
  
  // KernPtr k = *XPtr< KernPtr >( Kern );
  
  // Filling Gaussian process covariance matrix
  if ( symmetric ) {
    // #pragma omp parallel for private(i, j, x, y, Kern) shared( X, Y ) schedule( static ) collapse(2)
    // #pragma omp parallel for private( i, j, x, y ) shared( K )
    // #pragma omp parallel for shared( K )
#pragma omp parallel for shared( K )
    for ( i = 0; i < n; i++ ) { 
      for ( j = i; j < n; j++ ) {
        Eigen::VectorXd x, y;
        x = X.row( i );
        y = Y.row( j );
        // #pragma omp critical
        // {
        K( i, j ) = as<double>( Kern( x, y ) );
        // }
        
        if ( j > i ) {
          K( j, i ) = K( i, j );
        }
        
      }
    }
    
  } else { 
    // #pragma omp parallel for private(i, j, x, y) shared( X, Y ) schedule( static ) collapse(2)
    // #pragma omp parallel for private( i, j, x, y ) shared( K )
    // #pragma omp parallel for private( i, j ) shared( X, Y, K )
    for ( i = 0; i < n; i++ ) { 
      for ( j = 0; j < m; j++ ) {
        Eigen::VectorXd x, y;
        x = X.row( j );
        y = Y.row( i );
        
        // #pragma omp critical
        // {
        K( i, j ) = as<double>( Kern( x, y ) );
        // }
      }
    }
  }
  
  return K;
  
}

//--------------------------------------------------------------------------------------------------
List Krig( const Eigen::MatrixXd& Z, 
           const Eigen::MatrixXd& K, 
           const Eigen::MatrixXd& k,
           const Eigen::MatrixXd& G, 
           const Eigen::MatrixXd& g,
           const std::string type, 
           const std::string cinv ) {
  
  // The dimensions considerations
  // dim( K ) = n x n
  // dim( k_r ) = n x m
  // dim( l_r ) = dim( J * k_r )  = n x m
  // dim( Z ) = n x 1
  // dim( W ) = dim( k_r^t * J^t * Z ) = m x 1
  int n = Z.rows();
  int m = k.cols();
  
  List KRIG;
  
  Eigen::MatrixXd J( n, n );
  Eigen::MatrixXd W( m, n );
  Eigen::MatrixXd L( m, n );
  
  // Inverse computation
  if ( cinv == "syminv" ) {
    J = inv_sympd( K );
    
  } else if ( cinv == "inv" ) {
    J = inv( K );
    
  } else if ( cinv == "cholinv" ) {
    J = chol( K );
    J = inv( J );
    J = J.transpose() * J;
    
  } else if ( cinv == "ginv" ) {
    J = K;
  }
  
  // Kriging computation
  if ( type == "simple" ) { // Simple kriging
    
    L = J *  k;
    W = L.transpose() * Z;  
    
    KRIG[ "Z" ] = W;
    KRIG[ "L" ] = L;
    KRIG[ "J" ] = J;
    
    
  } else if ( type == "ordinary" ) {  // Ordinary kriging
    double alpha;
    Eigen::MatrixXd u = Eigen::MatrixXd::Ones( n, 1 );
    Eigen::MatrixXd tau( n, m );
    
    alpha = 1.0 / u.transpose() * J * u;
    tau = Eigen::MatrixXd::Ones( n, m ) - Eigen::MatrixXd::Ones( n, n ) * J * k;
    
    L = J * ( k + alpha * tau ) ;
    W = L.transpose() * Z;
    
    KRIG[ "Z" ] = W;
    KRIG[ "L" ] = L;
    KRIG[ "J" ] = J;
    KRIG[ "alpha" ] = alpha;
    KRIG[ "tau" ] = tau;
    
  } else if ( type == "universal" ) { // Universal kriging
    
    int p = G.rows();
    Eigen::MatrixXd A( n, p );
    Eigen::MatrixXd tau( p, m );
    
    A = G.transpose() * inv_sympd( G * J * G.transpose() );
    tau = g - G * J * k;
    
    L = J * ( k + A * tau ) ;
    W = L.transpose() * Z;
    
    KRIG[ "Z" ] = W;
    KRIG[ "L" ] = L;
    KRIG[ "J" ] = J;
    KRIG[ "A" ] = A;
    KRIG[ "tau" ] = tau;
    
  }
  
  return KRIG;
  
}

