#ifndef __KRIG_variog__
#define __KRIG_variog__

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <vector>


// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Computes the variogram.
//' @description This useful function is commonly employed in the study of isotropic stationary 
//' spatial processes. 
//' @param Z Vector of observations.
//' @param X Points matrix.
//' @param d Distance function.
//' @return Variogram vector.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' library( KRIG )
//' vignette( topic = 'copper_mining_2d', package = 'KRIG' )
//' @export
// [[Rcpp::export]]
List variogram( const arma::mat& Z,
                const arma::mat& X, 
                Function d );
#endif
