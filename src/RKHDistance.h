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
                         const arma::colvec& p );
#endif
