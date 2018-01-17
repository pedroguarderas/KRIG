#ifndef __KRIG_sens__
#define __KRIG_sens__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Combinatorial variance computation.
//' @description For a given combination this function computes the associated variance for the
//' variable enumerated by the combination values.
//' @param KF 
//' @param comb Combination.
//' @param X Points grid.
//' @param Gamma Cube with integral results.
//' @return Real value of sensitivity.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double sens_idx( const arma::colvec KF,
                 const arma::colvec comb,
                 const arma::mat X,
                 const arma::cube Gamma );


//--------------------------------------------------------------------------------------------------
//' @title Combinatorial variance computation.
//' @description Computation of variance
//' @param Gamma Cube with integral results.
//' @return Real value of sensitivity.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double sens_var( const arma::colvec KF, 
                 const arma::cube Gamma );

#endif
