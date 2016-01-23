/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 02-04-2013
  file: distance.hpp

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

//#ifndef DISTANCE
//#define DISTANCE

//#include <gsl/gsl_math.h>
#include <math.h>
#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;

/*__________________________________________________________________________________________________
  Distance functions
*/

//namespace krig {

// [[Rcpp::export]]
double dst_weighted( const NumericVector& x, const NumericVector& y, 
const NumericVector& w, const NumericVector& p ) {
  
  double d = 0.0;
  
  if ( x.size() > 0 && x.size() == y.size() && y.size() == w.size() && w.size() == p.size() ) {
    int i;
    #pragma omp parallel for shared( x, y, w, p ) private( i ) reduction (+:d)
    for( i = 0; i < x.size() ; i++ ) {
      d += w[i] * pow( abs( x[i] - y[i] ), p[i] );
    }
  }
  return d;
}

// [[Rcpp::export]]
double dst( const NumericVector& x, const NumericVector& y ) {
  
  double d = 0.0;
  
  if ( x.size() > 0 && x.size() == y.size() ) {
    int i;
    #pragma omp parallel for shared( x, y ) private( i ) reduction (+:d)
    for( i = 0; i < x.size() ; i++ ) {
      d += ( x[i] - y[i] ) * ( x[i] - y[i] );
    }
  }
  return d;
}


//} // namespace krig

//#endif // DISTANCE
