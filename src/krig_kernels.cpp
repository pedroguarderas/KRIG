
#include "krig_kernels.h"

//--------------------------------------------------------------------------------------------------
double linear_kernel( const arma::colvec& x, const arma::colvec& y, const double& alpha ) {
  return as_scalar( alpha * x.t() * y );
}

//--------------------------------------------------------------------------------------------------
double polynomial_kernel( const arma::colvec& x, const arma::colvec& y, 
                          const double& alpha, const double& beta, const double& n ) {
  return pow( as_scalar( alpha * x.t() * y ) + beta, n );
}

//--------------------------------------------------------------------------------------------------
double square_kernel( const double& h, const double& alpha ) {
  return alpha * h * h;
}

//--------------------------------------------------------------------------------------------------
double triangular_kernel( const double& h, const double& c, const double& alpha ) {
  return c * GSL_MAX_DBL( alpha - h, 0 );
}

//--------------------------------------------------------------------------------------------------
double exp_kernel( const double& h, const double& sigma, const double& theta ) {
  return sigma * gsl_sf_exp( -h / theta );
}

//--------------------------------------------------------------------------------------------------
double gaussian_kernel( const double& h, const double& sigma, const double& theta ) {
  return sigma * gsl_sf_exp( -h * h / theta );
}

//--------------------------------------------------------------------------------------------------
double spherical_kernel( const double& h, const double& phi, const double& theta ) {
  double Ker = 0.0;
  if ( h < theta ) {
    Ker = phi * ( 1 - 1.5 * h / theta + 0.5 * h * h * h / ( theta * theta * theta ) );
  } 
  return Ker;
}

//--------------------------------------------------------------------------------------------------
double matern_kernel( const double& h, const double& v, const double& sigma, const double& theta ) {
  double H = sqrt( v ) * h / theta;
  if ( H > 0 ) {
    return sigma * sigma * 2 * pow( H, v ) * gsl_sf_bessel_Knu( v, 2 * H ) / gsl_sf_gamma( v );
  }
  else {
    return 1.0;
  }
}

//--------------------------------------------------------------------------------------------------
double multilog_kernel( const double& h, const double& R ) {
  return gsl_sf_log( h * h + R * R );
}

//--------------------------------------------------------------------------------------------------
double nat_cubic_spline_kernel( const double& h, const double& R ) {
  return pow( h * h + R * R, 2.0 / 3.0 );
}

//--------------------------------------------------------------------------------------------------
double  thin_plate_kernel( const double& h, const double& R ) {
  return ( h * h + R * R ) * gsl_sf_log( h * h + R * R );
}
