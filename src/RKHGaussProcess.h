#ifndef __RKHGausProcess__
#define __RKHGausProcess__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title  Covariance kernel matrix
//' @description Computes covariance kernel matrix
//' @param X
//' @param Y
//' @param k Kernel
//' @param symmetric check if matrix will be symmetric
//' @return Covariance matrix
//' @author Pedro Guarderas
//' @useDynLib RKHSENS
//' @importFrom Rcpp sourceCpp
//' @exportPattern("^[[:alpha:]]+")
//' @export
// [[Rcpp::export]]
arma::mat RKHCov( const arma::mat& X, const arma::mat& Y, Function Kern, const bool symmetric = false );


//--------------------------------------------------------------------------------------------------
//' @title Gaussian regression
//' @description Computes Gaussian regression with a given covariance kernel
//' @param Z
//' @param X
//' @param Y 
//' @param Kern kernel
//' @return List
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List RKHGaussProcess( const arma::mat& Z, const arma::mat& X, const arma::mat& Y, Function Kern );

#endif
