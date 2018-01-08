#ifndef __RKHDistance__
#define __RKHDistance__

#include <RcppArmadillo.h>
#include <cmath>
#include <omp.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(oepnmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Generic weighted p-distance
//' @description Many isotropic models can be defined employing a distance. In particular this 
//' function implements a distance employing weights and different powers.
//' @param x First vector
//' @param y Second vector
//' @param w Weights for every coordinate in the vectors
//' @param p Powers for every coordinate in the vectors. To define a true norm every coordinate
//' has to be greater than 1.
//' @return Numeric value of the weighted p-distance.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHWeightPowDist( const arma::colvec& x, 
                         const arma::colvec& y, 
                         const arma::colvec& w, 
                         const arma::colvec& p );
#endif
