
#include "krig_sensitivity.h"

//--------------------------------------------------------------------------------------------------
double Krigidx( const Eigen::VectorXd& KF,
                const Eigen::VectorXd& comb,
                const Eigen::MatrixXd& X,
                const Eigen::ArrayXd& Gamma ) {
  int k;
  int n = X.n_rows;
  int c = comb.size();
  double S;
  
  Eigen::MatrixXd W = Eigen::MatrixXd::Ones( n, n );
  
  for ( k = 0; k < c; k++ ) {
    W = W % Gamma.slice( comb( k ) - 1 );
  }
  S = as_scalar( KF.t() * W * KF );
  
  return S;
}

//--------------------------------------------------------------------------------------------------
double Krigvar( const Eigen::VectorXd& KF, 
                const Eigen::ArrayXd& Gamma ) {
  
  int i;
  int n = Gamma.n_rows;
  int m = Gamma.n_slices;
  double Var;
  
  Eigen::MatrixXd V = Eigen::MatrixXd::Ones( n, n );
  for( i = 0; i < m; i++ ) {
    V = V % ( Eigen::MatrixXd::Ones( n, n ) + Gamma.slice( i ) );
  }
  V = V - Eigen::MatrixXd::Ones(( n, n );
  
  Var = as_scalar( KF.t() * V * KF );
  
  return Var;
}
