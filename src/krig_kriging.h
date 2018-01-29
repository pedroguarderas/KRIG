#ifndef __KRIG_kriging__
#define __KRIG_kriging__

#include <RcppArmadillo.h>
#include <string>
#include <omp.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(oepnmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//--------------------------------------------------------------------------------------------------
//' @title Spatial covariance matrix.
//' @description To compute a kriging, it is necessary the spatial covariance matrix. The spatial 
//'   covariance could computed between to sets of points X and Y with different dimension and the 
//'   result it is not necessarily a square matrix.
//' @param X First set of spatial points.
//' @param Y Second set of spatial points.
//' @param Kern Kernel function.
//' @param symmetric If result of computation will be a square matrix, the time computation can
//'   be improved setting this parameter to TRUE, the default is FALSE.
//' @return The spatial covariance matrix.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @useDynLib KRIG
//' @importFrom Rcpp sourceCpp
//' @exportPattern("^[[:alpha:]]+")
//' @seealso For a complete application you can check the documentation of \code{\link{Krig}}.
//' @export
// [[Rcpp::export]]
arma::mat Kov( const arma::mat& X, 
               const arma::mat& Y, 
               Function Kern, 
               const bool symmetric = false );

//--------------------------------------------------------------------------------------------------
//' @title Kriging computation.
//' @description Computes the kriging linear estimator for different types of kriging models.
//' @param Z Observed values of the spatial process.
//' @param K Covariance matrix computed for the position \eqn{X} where the spatial process \eqn{Z}
//' was observed.
//' @param k Covariance matrix computed for the position \eqn{X} where the spatial process \eqn{Z}
//' was observed and the position \eqn{Y} where the spatial process \eqn{Z} will be predicted.
//' @param G When universal kriging will be computed, this matrix represents the values of the 
//' of the functions representing the mean of the process \eqn{Z}, evaluated in the spatial 
//' points \eqn{X} where the spatial process was first observed.
//' @param g When universal kriging will be computed, this matrix represents the evaluation of the
//' functions representing the mean over the new position points \eqn{Y} where the spatial process
//' \eqn{Z} will be predicted. 
//' @param type Type of kriging model, possible values are: simple, ordinary, universal.
//' @param cinv Specifies how the inverse of the covariance matrix \eqn{K} will be computed. 
//' Possible values are: syminv = symmetric matrix inverse computation, inv = usual armadillo
//' inverse computation, cholinv = Cholesky based inverse computation, ginv = given inverse not
//' necessary to compute inverse at all.
//' @return Depending of the type of analysis the list of results change.
//' \item{Z}{New estimated values for Z.}
//' \item{L}{Linear coefficients determined by kriging.}
//' \item{J}{Inverse of the covariance matrix.}
//' \item{tau}{Factor computed in the ordinary and universal kriging.}
//' \item{alpha}{Factor computed in the ordinary kriging.}
//' \item{A}{Factor computed in the universal kriging.}
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' library( KRIG )
//' vignette( topic = 'simple_kriging', package = 'KRIG' )
//' vignette( topic = 'ordinary_kriging', package = 'KRIG' )
//' vignette( topic = 'universal_kriging', package = 'KRIG' )
//' vignette( topic = 'copper_mining_2d', package = 'KRIG' )
//' @export
// [[Rcpp::export]]
List Krig( const arma::mat& Z, 
           const arma::mat& K, 
           const arma::mat& k,
           const arma::mat& G, 
           const arma::mat& g,
           const std::string type = "ordinary", 
           const std::string cinv = "syminv" );

#endif
