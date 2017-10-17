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
double RKHKerLinear( const double& h, const double& alpha = 1.0 ) {
  return alpha * h;
}

//--------------------------------------------------------------------------------------------------
//' @title Square kernel
//' @description
//' @param h
//' @param alpha
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerSqr( const double& h, const double& alpha = 1.0 ) {
  return alpha * h * h;
}

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
double RKHKerTri( const double& h, const double& c = 1.0, const double& alpha = 1.0 ) {
  return c * GSL_MAX_DBL( alpha - h, 0 );
}

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
double RKHKerExp( const double& h, const double& sigma = 1.0, const double& theta = 1.0 ) {
  return sigma * sigma * gsl_sf_exp( -h / theta );
}

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
double RKHKerSqrExp( const double& h, const double& sigma = 1.0, const double& theta = 1.0 ) {
  return sigma * sigma * gsl_sf_exp( -h * h / ( theta * theta ) );
}

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
const double& theta = 1.0 ) {
  double H = sqrt( v ) * h / theta;
  if ( H > 0 ) {
    return sigma * sigma * 2 * pow( H, v ) * gsl_sf_bessel_Knu( v, 2 * H ) / gsl_sf_gamma( v );
  }
  else {
    return 1.0;
  }
}

//--------------------------------------------------------------------------------------------------
//' @title Multilog kernel
//' @description
//' @param h
//' @param R
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerMultilog( const double& h, const double& R = 1.0 ) {
  return gsl_sf_log( h * h + R * R );
}

//--------------------------------------------------------------------------------------------------
//' @title Natural cubic spline kernel
//' @description
//' @param h
//' @param R
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerNatCubSpl( const double& h, const double& R = 1.0 ) {
  return pow( h * h + R * R, 2.0 / 3.0 );
}

//--------------------------------------------------------------------------------------------------
//' @title Thin plate kernel
//' @description
//' @param h
//' @param R
//' @return Real value
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
double RKHKerPlateSpl( const double& h, const double& R = 1.0 ) {
  return ( h * h + R * R ) * gsl_sf_log( h * h + R * R );
}

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
double RKHKerMix( const double& h, const double& sigma = 1.0, const double& theta = 1.0 ) {
  double ht = h / theta;
  return sigma * sigma * ( 1 + ( sqrt(5.0) + 5.0 / 3.0 * ht ) * ht ) * 
  gsl_sf_exp( -sqrt(5.0) * ht );
}

#endif
