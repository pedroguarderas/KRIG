
#include "RKHKernels.h"

//--------------------------------------------------------------------------------------------------
double RKHKerLinear( const double& h, const double& alpha = 1.0 ) {
  return alpha * h;
}

//--------------------------------------------------------------------------------------------------
double RKHKerSqr( const double& h, const double& alpha = 1.0 ) {
  return alpha * h * h;
}

//--------------------------------------------------------------------------------------------------
double RKHKerTri( const double& h, const double& c = 1.0, const double& alpha = 1.0 ) {
  return c * GSL_MAX_DBL( alpha - h, 0 );
}

//--------------------------------------------------------------------------------------------------
double RKHKerExp( const double& h, const double& sigma = 1.0, const double& theta = 1.0 ) {
  return sigma * sigma * gsl_sf_exp( -h / theta );
}

//--------------------------------------------------------------------------------------------------
double RKHKerSqrExp( const double& h, const double& sigma = 1.0, const double& theta = 1.0 ) {
  return sigma * sigma * gsl_sf_exp( -h * h / ( theta * theta ) );
}

//--------------------------------------------------------------------------------------------------
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
double RKHKerMultilog( const double& h, const double& R = 1.0 ) {
  return gsl_sf_log( h * h + R * R );
}

//--------------------------------------------------------------------------------------------------
double RKHKerNatCubSpl( const double& h, const double& R = 1.0 ) {
  return pow( h * h + R * R, 2.0 / 3.0 );
}

//--------------------------------------------------------------------------------------------------
double RKHKerPlateSpl( const double& h, const double& R = 1.0 ) {
  return ( h * h + R * R ) * gsl_sf_log( h * h + R * R );
}

//--------------------------------------------------------------------------------------------------
double RKHKerMix( const double& h, const double& sigma = 1.0, const double& theta = 1.0 ) {
  double ht = h / theta;
  return sigma * sigma * ( 1 + ( sqrt(5.0) + 5.0 / 3.0 * ht ) * ht ) * 
    gsl_sf_exp( -sqrt(5.0) * ht );
}
