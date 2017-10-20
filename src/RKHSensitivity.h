#ifndef __RKHCombination__
#define __RKHCombination__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Sensitiviy analysis
//' @description Computation of Sobol index
//' @param KF 
//' @param comb Combination.
//' @param X Points grid.
//' @param Gamma Cube with integral results.
//' @return Real value of sensitivity.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHSobolIndex( const arma::colvec KF,
                      const arma::colvec comb,
                      const arma::mat X,
                      const arma::cube Gamma );


//--------------------------------------------------------------------------------------------------
//' @title Var
//' @description Computation of variance
//' @param Gamma Cube with integral results.
//' @return Real value of sensitivity.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHSobolVar( const arma::colvec KF, 
                    const arma::cube Gamma );

#endif
