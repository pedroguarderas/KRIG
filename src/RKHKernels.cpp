
#include "RKHKernels.h"

//--------------------------------------------------------------------------------------------------
double RKHKerLinear( const double& h, const double& alpha ) {
  return alpha * h;
}

//--------------------------------------------------------------------------------------------------
double RKHKerSqr( const double& h, const double& alpha ) {
  return alpha * h * h;
}

//--------------------------------------------------------------------------------------------------
double RKHKerTri( const double& h, const double& c, const double& alpha ) {
  return c * GSL_MAX_DBL( alpha - h, 0 );
}

//--------------------------------------------------------------------------------------------------
double RKHKerExp( const double& h, const double& sigma, const double& theta ) {
  return sigma * sigma * gsl_sf_exp( -h / theta );
}

//--------------------------------------------------------------------------------------------------
double RKHKerSqrExp( const double& h, const double& sigma, const double& theta ) {
  return sigma * sigma * gsl_sf_exp( -h * h / ( theta * theta ) );
}

//--------------------------------------------------------------------------------------------------
double RKHKerSpher( const double& h, const double& phi, const double& theta ) {
  double Ker = 0.0;
  if ( h < theta ) {
    Ker = phi * ( 1 - 1.5 * h / theta + 0.5 * h * h * h / ( theta * theta * theta ) );
  } 
  return Ker;
}

//--------------------------------------------------------------------------------------------------
double RKHKerMatern( const double& h, const double& v, const double& sigma, const double& theta ) {
  double H = sqrt( v ) * h / theta;
  if ( H > 0 ) {
    return sigma * sigma * 2 * pow( H, v ) * gsl_sf_bessel_Knu( v, 2 * H ) / gsl_sf_gamma( v );
  }
  else {
    return 1.0;
  }
}

//--------------------------------------------------------------------------------------------------
double RKHKerMultilog( const double& h, const double& R ) {
  return gsl_sf_log( h * h + R * R );
}

//--------------------------------------------------------------------------------------------------
double RKHKerNatCubSpl( const double& h, const double& R ) {
  return pow( h * h + R * R, 2.0 / 3.0 );
}

//--------------------------------------------------------------------------------------------------
double RKHKerPlateSpl( const double& h, const double& R ) {
  return ( h * h + R * R ) * gsl_sf_log( h * h + R * R );
}

//--------------------------------------------------------------------------------------------------
double RKHKerMix( const double& h, const double& sigma, const double& theta ) {
  double ht = h / theta;
  return sigma * sigma * ( 1 + ( sqrt(5.0) + 5.0 / 3.0 * ht ) * ht ) * 
    gsl_sf_exp( -sqrt(5.0) * ht );
}
