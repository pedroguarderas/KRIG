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

/*__________________________________________________________________________________________________
  Numerical integration funtion
 */

// [[Rcpp::export]]
double integralRKH( NumericVector x, Function k ) {
  double I = 0.0;
  int i;
  
  #pragma omp parallel for private( i ) reduction (+:I)
  for( i = 0; i < x.size() - 1; i++ ) {
    I += as<double>( k( x[i] ) ) * ( x[i+1] - x[i] );
  }
  return I;
}

