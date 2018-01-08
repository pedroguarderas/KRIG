#ifndef __RKHEstimate__
#define __RKHEstimate__

#include <RcppArmadillo.h>
#include <string>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Spatial covariance matrix.
//' @description To compute a Gaussian regression or interpolation is necessary use Kernel, and
//'   it is requeried to construct the spatial convariance matrix. The spatial covariance could
//'   be estimated between to sets of points X and Y and the result is not necessarily a square
//'   matrix.
//' @param X First set of spatial points.
//' @param Y Second set of spatial points.
//' @param Kern Kernel function.
//' @param symmetric If result of computation will be a square matrix, the time computation can
//'   be improved setting this parameter to TRUE, the default is FALSE.
//' @return The spatial covariance matrix.
//' @author Pedro Guarderas
//' @useDynLib RKHSENS
//' @importFrom Rcpp sourceCpp
//' @exportPattern("^[[:alpha:]]+")
//' @export
// [[Rcpp::export]]
arma::mat RKHCov( const arma::mat& X, 
                  const arma::mat& Y, 
                  Function Kern, 
                  const bool symmetric = false );

//--------------------------------------------------------------------------------------------------
//' @title Gaussian regression
//' @description Computes Gaussian regression with a given covariance kernel
//' @param Z Observed spatial values.
//' @param X Spatial points, ubication where the Z value was observed.
//' @param Y Spatial points where it is requested to predict Z.
//' @param Kern Kernel functions.
//' @param type Type of Gaussian regression, 0 = simple interpolation, 1 = simple regression, 
//'  2 = ordinary interpolation, 3 = ordinary regression, 4 = universal interpolation, 
//'  5 = universal regression.
//' @param cinv Form to compute the inverse.
//' @return List containing the different important matrices employed for the Gaussian regression.
//' \item{W}{Predicted new values for Z}
//' \item{K}{Covariance matrix for X}
//' \item{k}{Covariance matrix between Y and X}
//' \item{J}{Inverse of the covariance matrix}
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List RKHEstimate( const arma::mat& Z, 
                  const arma::mat& K, 
                  const arma::mat& k,
                  const arma::mat& G, 
                  const arma::mat& g,
                  const std::string type = "ordinary", 
                  const std::string typeinv = "syminv" );

#endif
