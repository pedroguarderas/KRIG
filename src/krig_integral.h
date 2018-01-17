#ifndef __KRIG_integral__
#define __KRIG_integral__

#include <RcppArmadillo.h>
#include "krig_kriging.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title One coordinate kernel integral.
//' @description This function is part of the routines employed in the sensitivity analysis, takes 
//' the kernel \eqn{k} and for each fixed coordinate in \eqn{x}, the integral in the second 
//' variable \eqn{y}, is computed in the interval \eqn{a} to \eqn{b} by taking \eqn{n} uniform
//' steps.
//' @param Kern Kernel function.
//' @param x Column vector with values for the first coordinate of the kernel.
//' @param a Inferior limit for the integral in y.
//' @param b Superior limit for the integral in y.
//' @param n Number of uniform division to compute the integral.
//' @return Vector with integrals while the x coordinate is fixed.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @export
// [[Rcpp::export]]
arma::colvec vector_integrate_kernel( Function Kern, const arma::colvec x, 
                                      const double& a, const double& b, const double& n );

//--------------------------------------------------------------------------------------------------
//' @title Complete kernel integral.
//' @description This function is part of the routines employed in the sensitivity analysis, 
//' calculate the integral in both coordinate \eqn{x} and \eqn{y} of the kernel, over the square
//' domain give by the limits \eqn{a} and \eqn{b}. 
//' @param Kern Kernel function.
//' @param a Inferior limit for the integral in each coordinate.
//' @param b Superior limit for the integral in each coordinate.
//' @param n Number of uniform division to compute the integral.
//' @return Real value with the integral value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @export
// [[Rcpp::export]]
double complete_integrate_kernel( Function Kern, const double& a, const double& b, const double& n );

//--------------------------------------------------------------------------------------------------
//' @title Integrals of a list of kernels.
//' @description This function is part of the routines employed in the sensitivity analysis, 
//' computes vector of integrals and complete integrals of kernels specified in the data frame
//' Kernels.
//' @param Kernels data.frame of kernels composed by four columns.
//' \enumerate{
//'   \item Kernel name.
//'   \item Inferior limit for integral.
//'   \item Superior limit for integral.
//'   \item Number of steps for discretization of integrals.
//' }
//' @param X matrix containing in each row the coordinate where the one coordinate integrals will 
//' be evaluated.
//' @return List with one coordinate integrals and complete kernel integrals.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @export
// [[Rcpp::export]]
List list_integrate_kernel( const DataFrame& Kernels, const arma::mat& X );

//--------------------------------------------------------------------------------------------------
//' @title KANOVA, kernel anova under RKHS approximations.
//' @description Under an approximation to the sensitivity analysis based in variance computation
//' the different indexes of combinatorial sensitivity values can be computed employing the
//' values of kernel integrals.
//' @param Kernels data.frame of kernels composed by four columns.
//' \enumerate{
//'   \item Kernel name.
//'   \item Inferior limit for integral.
//'   \item Superior limit for integral.
//'   \item Number of steps for discretization of integrals.
//' }
//' @param Integral A list containing the results of kernel integrals of the functions
//' \code{\link{vector_integrate_kernel}}.
//' @param X matrix containing in each row the coordinate where the one coordinate integrals will 
//' be evaluated.
//' @return List with contaning the Gamma 3D array where the different combination variance are
//' stocked and the total matrix variance named Kanova. 
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @export
// [[Rcpp::export]]
List Kanova( const DataFrame& Kernels, const List& Integral,  const arma::mat& X );

#endif
