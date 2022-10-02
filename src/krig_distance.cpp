
#include "krig_distance.h"

double weight_pow_dist( const Eigen::VectorXd& x, 
                        const Eigen::VectorXd& y, 
                        const Eigen::VectorXd& w, 
                        const Eigen::VectorXd& p ) {
  
  double d = 0.0;
  
  if ( x.size() > 0 && x.size() == y.size() && y.size() == w.size() && w.size() == p.size() ) {
    int i, n;
    n = x.size();
    #pragma omp parallel for shared( w, x, y, p ) private( i ) reduction(+:d)
    for( i = 0; i < n; i++ ) {
      d += w(i) * std::pow( std::abs( x(i) - y(i) ), p(i) );
    }
  }
  return d;
}
