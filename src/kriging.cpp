/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 12-04-2015
  file: kriging.cpp

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
List krigingSimpleRKH( NumericMatrix Z, NumericMatrix X, NumericMatrix x0, Function k ) {
  int i, j, l;
  int n = X.nrow();
  int m = x0.nrow();
  List Kriging;
  NumericMatrix k0( m, n );
  NumericMatrix K( n, n );
  
  for ( i = 0; i < n; i++ ) { 
    for ( j = i; j < n; j++ ) {
      
      K(i,j) = as<double>( k( X( i, _ ), X( j, _ ) ) );
      
      if ( j > i ) {
        K(j,i) = K(i,j);
      }
    }
    
    for ( l = 0; l < m; l++ ) {
      k0(l,i) = as<double>( k( x0( l, _ ), X( i, _ ) ) );
    }
  }
  
//  L<-chol( K )
//  J<-chol2inv( L )
//  Z0<-k0 %*% J %*% Z
  

//  Kriging["Z0"] = Z0;
  Kriging["K"] = K;
  Kriging["k0"] = k0;
//  Kriging["L"] = L;
//  Kriging["J"] = J;

  return( Kriging );
}

