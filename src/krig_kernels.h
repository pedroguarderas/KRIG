#ifndef __KRIG_kernels__
#define __KRIG_kernels__

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
 Anisotropic.
 In the current section are defined different anisotropic kernels.
 ---------------------------------------------------------------------------------------------------
*/

//--------------------------------------------------------------------------------------------------
//' @title Linear kernel
//' @description Anisotropic kernel defined by the scalar product.
//' @param x first column vector.
//' @param y second column vector.
//' @param alpha amplitude parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' n<-10
//' x<-matrix( runif( n ), n, 1 )
//' y<-matrix( runif( n ), n, 1 )
//' alpha<-0.5
//' linear_kernel( x, y, alpha )
//' @export
// [[Rcpp::export]]
double linear_kernel( const arma::colvec& x, const arma::colvec& y, const double& alpha );

//--------------------------------------------------------------------------------------------------
//' @title Polynomial kernel
//' @description Anisotropic kernel defined like a polynomial in the scalar product.
//' @param x first column vector.
//' @param y second column vector.
//' @param alpha amplitude parameter.
//' @param beta displacement parameter.
//' @param n power parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' N<-10
//' x<-matrix( runif( N ), N, 1 )
//' y<-matrix( runif( N ), N, 1 )
//' alpha<-0.5
//' beta<-2
//' n<-3
//' polynomial_kernel( x, y, alpha, beta, n )
//' @export
// [[Rcpp::export]]
double polynomial_kernel( const arma::colvec& x, const arma::colvec& y, 
                          const double& alpha, const double& beta, const double& n );


/*--------------------------------------------------------------------------------------------------
 Isotropic.
 In the current section are defined different isotropic kernels.
 ---------------------------------------------------------------------------------------------------
*/

//--------------------------------------------------------------------------------------------------
//' @title Square kernel.
//' @description Isotropic kernel given by the square distance.
//' @param h distance variable.
//' @param alpha amplitude parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-4
//' alpha<-2
//' square_kernel( h, alpha )
//' @export
// [[Rcpp::export]]
double square_kernel( const double& h, const double& alpha = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Triangular kernel.
//' @description Isotropic kernel defined with the max function.
//' @param h distance variable.
//' @param c amplitude parameter.
//' @param alpha maximum distance value.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-0.2
//' c<-5
//' alpha<-2
//' triangular_kernel( h, c, alpha )
//' @export
// [[Rcpp::export]]
double triangular_kernel( const double& h, const double& c = 1.0, const double& alpha = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Exponential kernel.
//' @description Isotropic kernel.
//' @param h distance variable.
//' @param sigma amplitude parameter.
//' @param theta smoothness parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-0.3
//' sigma<-1.0
//' theta<-10
//' exp_kernel( h, sigma, theta )
//' @export
// [[Rcpp::export]]
double exp_kernel( const double& h, const double& sigma = 1.0, const double& theta = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Gaussian kernel.
//' @description Isotropic kernel.
//' @param h distance variable.
//' @param sigma amplitude parameter.
//' @param theta smoothness parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-4
//' alpha<-2
//' gaussian_kernel( h, alpha )
//' @export
// [[Rcpp::export]]
double gaussian_kernel( const double& h, const double& sigma = 1.0, const double& theta = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Spherical kernel.
//' @description Isotropic kernel commonly employed in geostatistic.
//' @param h distance variable.
//' @param phi amplitude parameter.
//' @param theta smoothness parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-4.0
//' phi<-1.0
//' theta<-100.0
//' spherical_kernel( h, phi, theta )
//' @export
// [[Rcpp::export]]
double spherical_kernel( const double& h, const double& phi, const double& theta );

//--------------------------------------------------------------------------------------------------
//' @title Matérn kernel.
//' @description Isotropic kernel.
//' @param h distance variable.
//' @param v power parameter.
//' @param sigma amplitude parameter.
//' @param theta smoothness parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-4.0
//' v<-2.0
//' sigma<-2.0
//' theta<-100.0
//' matern_kernel( h, v, sigma, theta )
//' @export
// [[Rcpp::export]]
double matern_kernel( const double& h, const double& v = 2.0, const double& sigma = 1.0,
                      const double& theta = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Multilog kernel.
//' @description Isotropic kernel.
//' @param h distance parameter.
//' @param R displacement parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-4
//' R<-2.3
//' multilog_kernel( h, R )
//' @export
// [[Rcpp::export]]
double multilog_kernel( const double& h, const double& R = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Natural cubic spline kernel.
//' @description Isotropic kernel.
//' @param h distance variable.
//' @param R displacement parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-4.0
//' R<-2.3
//' nat_cubic_spline_kernel( h, R )
//' @export
// [[Rcpp::export]]
double nat_cubic_spline_kernel( const double& h, const double& R = 1.0 );

//--------------------------------------------------------------------------------------------------
//' @title Thin plate kernel.
//' @description Isotropic kernel.
//' @param h distance variable. 
//' @param R displacement parameter.
//' @return Real value.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' h<-4.0
//' R<-0.5
//' thin_plate_kernel( h, R )
//' @export
// [[Rcpp::export]]
double thin_plate_kernel( const double& h, const double& R = 1.0 );

#endif
