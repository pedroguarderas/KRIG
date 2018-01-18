
#include "krig_variogram.h"

//--------------------------------------------------------------------------------------------------
List variogram( const arma::mat& Z,
                const arma::mat& X, 
                Function d,
                const double& delta ) {
  int i, j, k;
  int n = X.n_rows;
  int N = n * ( n + 1 ) / 2;
  
  std::vector< int > I( N );
  arma::rowvec x, y;
  arma::colvec V( N );
  arma::colvec D( N );
  arma::colvec S( N );
  
  List Varg;
  
  for ( i = 0; i < N; i++ ) { 
    I[i] = i;
  }
  
  k = 0;
  for ( i = 0; i < n; i++ ) { 
    for ( j = i; j < n; j++ ) {
      x = X.row( i );
      y = X.row( j );
      
      D(k) = as<double>( d( x, y ) );
      S(k) = ( Z(i,0) - Z(j,0) ) * ( Z(i,0) - Z(j,0) );
    
      k++;
    }
  }
  
  sort( I.begin(), I.end(),
       [&]( const int& k, const int& l ) {
         return ( D(k) < D(l) );
       }
  );
  
  V( 0 ) = S( I[0] );
  for ( k = 1; k < N; k++ ) {
    V( k ) = ( V( k - 1 ) * k + S( I[ k ] ) ) / ( k + 1.0 );
  }
  
  Varg[ "variogram" ] = V;
  Varg[ "distance" ] = D;
  Varg[ "sort" ] = I;
  
  return Varg;
}
