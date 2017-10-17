#ifndef __RKHGausProcess__
#define __RKHGausProcess__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title  Covariance kernel matrix
//' @description Computes covariance kernel matrix
//' @param X
//' @param Y
//' @param k Kernel
//' @param symmetric check if matrix will be symmetric
//' @return Covariance matrix
//' @author Pedro Guarderas
//' @useDynLib RKHSENS
//' @importFrom Rcpp sourceCpp evalCpp
//' @exportPattern("^[[:alpha:]]+")
//' @export
// [[Rcpp::export]]
arma::mat RKHCov( const arma::mat& X, const arma::mat& Y, Function Kern, const bool symmetric = false ) {
  int i, j;
  int m = X.n_rows;
  int n = Y.n_rows;
  arma::mat K( m, n );
  
  // Filling Gaussian process covariance matrix
  if ( symmetric ) {
    for ( i = 0; i < m; i++ ) { 
      for ( j = i; j < n; j++ ) {
        
        K( i, j ) = as<double>( Kern( X.row( i ), Y.row( j ) ) );
        
        if ( j > i ) {
          K( j, i ) = K( i, j );
        }
        
      }
    }
    
  } else { 
    for ( i = 0; i < m; i++ ) { 
      for ( j = 0; j < n; j++ ) {
        K( i, j ) = as<double>( Kern( X.row( i ), Y.row( j ) ) );
      }
    }
  }
  
  return K;

}

//--------------------------------------------------------------------------------------------------
//' @title Gaussian regression
//' @description Computes Gaussian regression with a given covariance kernel
//' @param Z
//' @param X
//' @param Y 
//' @param Kern kernel
//' @return List
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List RKHGaussProcess( const arma::mat& Z, const arma::mat& X, const arma::mat& Y, Function Kern ) {

  int n = X.n_rows;
  int m = Y.n_rows;
  int p = Z.n_cols;
  
  arma::mat k( m, n );
  arma::mat K( n, n );
  arma::mat J( n, n );
  arma::mat W( m, p );
  arma::mat U( n, p );
  arma::mat O = arma::ones( n, 1 );
  
  K = RKHCov( X, X, Kern, true );
  k = RKHCov( Y, X, Kern );
  
  J = inv_sympd( K );
  // W = k * J * Z;
  
  U = ( 1.0 / ( O.t() * J * O ) ) * ( O.t() * J * Z );
  W = U + k * J * ( Z - U );

  return List::create( Named( "W" ) = W,
                       Named( "K" ) = K,
                       Named( "k" ) = k,
                       Named( "U" ) = U,
                       Named( "J" ) = J );
}

#endif
