#ifndef __KRIG_cokriging__
#define __KRIG_cokriging__

#include <RcppEigen.h>
#include <string>
#include <omp.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


//--------------------------------------------------------------------------------------------------
//' @title co-Kriging computation.
//' @description Computes the co-kriging linear estimator for different types models.
//' @param Z Matrix of observed values of the spatial process.
//' @param K Covariance matrix computed for the position \eqn{X} where the spatial process \eqn{Z}
//' was observed.
//' @param k Covariance cube computed for the position \eqn{X} where the spatial process \eqn{Z}
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
//' @export
// [[Rcpp::export]]
List coKrig( const Eigen::MatrixXd& Z, 
             const Eigen::MatrixXd& K, 
             const Eigen::ArrayXd& k,
             const Eigen::MatrixXd& G, 
             const Eigen::ArrayXd& g,
             const std::string type = "ordinary", 
             const std::string cinv = "syminv" );

#endif
