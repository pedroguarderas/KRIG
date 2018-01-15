
#include "krig_distance.h"

double weight_pow_dist( const arma::colvec& x, 
                        const arma::colvec& y, 
                        const arma::colvec& w, 
                        const arma::colvec& p ) {
  
  double d = 0.0;
  
  if ( x.size() > 0 && x.size() == y.size() && y.size() == w.size() && w.size() == p.size() ) {
    int i;
    #pragma omp parallel for shared( w, x, y, p ) private( i ) reduction(+:d)
    for( i = 0; i < x.size() ; i++ ) {
      d += w(i) * std::pow( std::abs( x(i) - y(i) ), p(i) );
    }
  }
  return d;
}
