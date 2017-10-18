#ifndef __RKHIntegral__
#define __RKHIntegral__

#include <RcppArmadillo.h>
#include "RKHEstimate.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Vector integral
//' @description Computes covariance kernel matrix
//' @param Kern
//' @param x
//' @param a
//' @param b
//' @param n
//' @return Vector with integrals
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
arma::colvec RKHIntegrateKern( Function Kern, const arma::colvec x, 
                               const double& a, const double& b, const double& n );

//--------------------------------------------------------------------------------------------------
//' @title Complete kernel integral
//' @description Complete kernel integral
//' @param Kern
//' @param a
//' @param b
//' @param n
//' @return Real
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHCompIntegKern( Function Kern, const double& a, const double& b, const double& n );

//--------------------------------------------------------------------------------------------------
//' @title Integrals of kernels
//' @description Compute integrals of kernels
//' @param Kernels List of kernels
//' @return List with vector integrals and complete kernel integrals
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List RKHKernInteg( const DataFrame& Kernels, const arma::mat& X );

//--------------------------------------------------------------------------------------------------
//' @title Integrals of kernels
//' @description Compute integrals of kernels
//' @param Kernels List of kernels
//' @return List with vector integrals and complete kernel integrals
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List RKHAnova( const DataFrame& Kernels, const List& Integral,  const arma::mat& X );

#endif
