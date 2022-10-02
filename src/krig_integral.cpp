
#include "krig_integral.h"


//--------------------------------------------------------------------------------------------------
Eigen::VectorXd vector_integrate_kernel( Function Kern, 
                                      const Eigen::VectorXd x, 
                                      const double& a, 
                                      const double& b, 
                                      const double& n ) {
  int i, k;
  int m = x.size();
  double h = ( b - a ) / ( n - 1 );
  Eigen::VectorXd I = Eigen::zeros( m );
  
  for ( k = 0; k < m; k++ ) {
    for ( i = 0; i < n; i++ ) {
      I( k ) += as< double >( Kern( x( k ), a + i * h ) ) * h;
    }
  }
  
  return I;
}

//--------------------------------------------------------------------------------------------------
double integrate_kernel( Function Kern, const double& a, const double& b, const double& n ) {
  int i, j;
  double h = ( b - a ) / ( n - 1 );
  double I = 0;
  
  for ( j = 0; j < n; j++ ) {
    for ( i = 0; i < n; i++ ) {
      I += as< double >( Kern( a + i * h, a + j * h ) ) * h * h;
    }
  }
  
  return I;
}

//--------------------------------------------------------------------------------------------------
List list_integrate_kernel( const DataFrame& Kernels, 
                            const Eigen::MatrixXd& X ) {
  
  int k;
  int N = Kernels.nrows();
  Eigen::MatrixXd integral( X.n_rows, N );
  Eigen::VectorXd alpha( N );
  
  Eigen::VectorXd a = Kernels[1];
  Eigen::VectorXd b = Kernels[2];
  Eigen::VectorXd n = Kernels[3];
  CharacterVector FName = Kernels[0];
  
  for ( k = 0; k < N; k++ ) {
    Function Kern( as< std::string >( FName[k] ) ); 
    
    integral.col( k ) = vector_integrate_kernel( Kern, X.col( k ), a( k ), b( k ), n( k ) );
    alpha( k ) = integrate_kernel( Kern, a( k ), b( k ), n( k ) );
    
  }
  
  return List::create( Named( "integral" ) = integral, 
                       Named( "alpha" ) = alpha );
}

//--------------------------------------------------------------------------------------------------
List Kanova( const DataFrame& Kernels, 
             const List& Integral, 
             const Eigen::MatrixXd& X ) {
  
  int N, n, i, j, k;
  
  N = Kernels.nrows();
  n = X.n_rows;
  CharacterVector FName = Kernels[0];

  Eigen::MatrixXd Kanova = Eigen::MatrixXd::Ones( n, n );
  Eigen::ArrayXd Gamma( n, n, N );
  Eigen::MatrixXd I = Integral[0];
  Eigen::VectorXd a = Integral[1];
  
  for ( k = 0; k < N; k++ ) {
    
    Function Kern( as< std::string >( FName[k] ) ); 

    for ( i = 0; i < n; i++ ) {
      for ( j = i; j < n; j++ ) {
        
        Gamma( i, j, k ) = as<double>( Kern( X( i, k ), X( j, k ) ) ) - I( i, k ) * I( j, k ) / a( k );
        Kanova( i, j ) = Kanova( i, j ) * ( 1.0 + Gamma( i, j, k ) );
        
        if ( j > i ) {
          Gamma( j, i, k ) = Gamma( i, j, k );
          Kanova( j, i ) = Kanova( i, j );
        }
      }
    }

  }
  return List::create( Named( "Gamma" ) = Gamma, 
                       Named( "Kanova" ) = Kanova );
}
