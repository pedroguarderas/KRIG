
#include "RKHDistance.h"

double RKHWeightPowDist( const arma::colvec& x, 
                         const arma::colvec& y, 
                         const arma::colvec& w, 
                         const arma::colvec& p ) {
  
  double d = 0.0;
  
  if ( x.size() > 0 && x.size() == y.size() && y.size() == w.size() && w.size() == p.size() ) {
    int i;
    for( i = 0; i < x.size() ; i++ ) {
      d += w(i) * std::pow( std::abs( x(i) - y(i) ), p(i) );
    }
  }
  return d;
}
