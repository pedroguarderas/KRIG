
#include "krig_sensitivity.h"

//--------------------------------------------------------------------------------------------------
double sens_idx( const arma::colvec KF,
                 const arma::colvec comb,
                 const arma::mat X,
                 const arma::cube Gamma ) {
  int k;
  int n = X.n_rows;
  int c = comb.size();
  double S;
  
  arma::mat W = arma::ones( n, n );
    
  for ( k = 0; k < c; k++ ) {
    W = W % Gamma.slice( comb( k ) - 1 );
  }
  S = as_scalar( KF.t() * W * KF );

  return S;
}

//--------------------------------------------------------------------------------------------------
double sens_var( const arma::colvec KF, 
                 const arma::cube Gamma ) {
  
  int i;
  int n = Gamma.n_rows;
  int m = Gamma.n_slices;
  double Var;
  
  arma::mat V = arma::ones( n, n );
  for( i = 0; i < m; i++ ) {
    V = V % ( arma::ones( n, n ) + Gamma.slice( i ) );
  }
  V = V - arma::ones( n, n );
  
  Var = as_scalar( KF.t() * V * KF );
  
  return Var;
}
