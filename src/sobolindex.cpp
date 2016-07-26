/*__________________________________________________________________________________________________
 
 autor: Pedro Guarderas
 email: ajusworkopensource@gmail.com
 date: 25-07-2016
 file: sobolindex.cpp
 
 This program is free software; you can redistribute it and/or modify it under the 
 terms of the GNU General Public License as published by the Free Software Foundation; 
 either version 2 of the License, or (at your option) any later version.
 __________________________________________________________________________________________________
 */

#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;
using namespace std;

/*__________________________________________________________________________________________________
  Sobol index
 */
// [[Rcpp::export]]
List SobolIndex( List Kernels, List GK, NumericVector F, NumericVector C ) {
  
  int I, N, M, i, j, k, l;
  double Var;
  List S, kernel, G;
  NumericVector X, Int;
  
  I = Kernels.size();
  kernel = Kernels[0];
  X = kernel["grid"];
  N = X.size();
  
  NumericMatrix CMB;
  NumericMatrix U( N, N ), V( N, N ), W( N, N ), K( N, N );
  
  fill( U.begin(), U.end(), 1.0 );
  
  K = solve( GK["KANOVA"], F );
  
  G = GK["Gamma"];
  V = U;
  for( i = 0; i < I; i++ ) {
    V = V * ( U + G[i] );
  }
  V = V - U;
  
  i = 0;
  for ( j = 0; j < C.size(); j++ ) {
    CMB = combn( I, C[j] );
    for ( k = 0; k < M; k++ ) {
      W = U;
      for ( l = 0; l < M; l++ ) {
        W = W * G[ CMB(l,k) - 1 ];
      }
      S[i] = t(K) %*% W %*% K;
      i++;
      // names(S)[length(S)]<-paste( 'C.', paste( CMB[,k], collapse='.' ), sep = '' )
    }
  }
  Var = (K) %*% V %*% K;
  return List::create( Named( "S" ) = S, Named( "Var" ) = Var );
}