#ifndef __RKHKernels__
#define __RKHKernels__

#include <RcppArmadillo.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/*--------------------------------------------------------------------------------------------------
  Isotropic covariance models.
  In the current section are defined different isotropic kernels.
--------------------------------------------------------------------------------------------------*/


//--------------------------------------------------------------------------------------------------
//' @title Linear kernel
//' @description
//' @param h
//' @param alpha
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerLinear( const double& h, const double& alpha = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Square kernel
//' @description
//' @param h
//' @param alpha
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerSqr( const double& h, const double& alpha = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Triangular kernel
//' @description
//' @param h
//' @param c
//' @param alpha
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerTri( const double& h, const double& c = 1.0, const double& alpha = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Exponential kernel
//' @description
//' @param h
//' @param sigma
//' @param theta
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerExp( const double& h, const double& sigma = 1.0, const double& theta = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Gaussian kernel
//' @description
//' @param h
//' @param sigma
//' @param theta
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerSqrExp( const double& h, const double& sigma = 1.0, const double& theta = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Matérn kernel
//' @description
//' @param h
//' @param v
//' @param sigma
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerMatern( const double& h, const double& v = 2.0, const double& sigma = 1.0,
const double& theta = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Multilog kernel
//' @description
//' @param h
//' @param R
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerMultilog( const double& h, const double& R = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Natural cubic spline kernel
//' @description
//' @param h
//' @param R
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerNatCubSpl( const double& h, const double& R = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Thin plate kernel
//' @description
//' @param h
//' @param R
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerPlateSpl( const double& h, const double& R = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Mix kernel
//' @description
//' @param h
//' @param sigma
//' @param theta
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerMix( const double& h, const double& sigma = 1.0, const double& theta = 1.0 );

#endif
