#ifndef __RKHDistance__
#define __RKHDistance__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Generic function to compute distances
//' @description This function computes a weighted distance between vectors x and y
//' @param x
//' @param y
//' @param w weights
//' @param p power terms
//' @return Numeric value with the distance
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHWeightPowDist( const arma::colvec& x, 
                         const arma::colvec& y, 
                         const arma::colvec& w, 
                         const arma::colvec& p ) {
  
  double d = 0.0;
  
  if ( x.size() > 0 && x.size() == y.size() && y.size() == w.size() && w.size() == p.size() ) {
    int i;
    for( i = 0; i < x.size() ; i++ ) {
      d += w(i) * pow( abs( x(i) - y(i) ), p(i) );
    }
  }
  return d;
}

#endif
