
#include "krig_integral.h"


//--------------------------------------------------------------------------------------------------
arma::colvec vector_integrate_kernel( Function Kern, const arma::colvec x, 
                                      const double& a, const double& b, const double& n ) {
  int i, k;
  int m = x.size();
  double h = ( b - a ) / ( n - 1 );
  arma::colvec I = arma::zeros( m );
  
  for ( k = 0; k < m; k++ ) {
    for ( i = 0; i < n; i++ ) {
      I( k ) += as< double >( Kern( x( k ), a + i * h ) ) * h;
    }
  }
  
  return I;
}

//--------------------------------------------------------------------------------------------------
double complete_integrate_kernel( Function Kern, const double& a, const double& b, const double& n ) {
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
List list_integrate_kernel( const DataFrame& Kernels, const arma::mat& X ) {
  
  int k;
  int N = Kernels.nrows();
  arma::mat integral( X.n_rows, N );
  arma::colvec alpha( N );
  
  arma::colvec a = Kernels[1];
  arma::colvec b = Kernels[2];
  arma::colvec n = Kernels[3];
  CharacterVector FName = Kernels[0];
  
  for ( k = 0; k < N; k++ ) {
    Function Kern( as< std::string >( FName[k] ) ); 
    
    integral.col( k ) = vector_integrate_kernel( Kern, X.col( k ), a( k ), b( k ), n( k ) );
    alpha( k ) = complete_integrate_kernel( Kern, a( k ), b( k ), n( k ) );
    
  }
  
  return List::create( Named( "integral" ) = integral, 
                       Named( "alpha" ) = alpha );
}

//--------------------------------------------------------------------------------------------------
List Kanova( const DataFrame& Kernels, const List& Integral, const arma::mat& X ) {
  
  int N, n, i, j, k;
  
  N = Kernels.nrows();
  n = X.n_rows;
  CharacterVector FName = Kernels[0];

  arma::mat Kanova = arma::ones( n, n );
  arma::cube Gamma( n, n, N );
  arma::mat I = Integral[0];
  arma::colvec a = Integral[1];
  
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
