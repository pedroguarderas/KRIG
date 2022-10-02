#ifndef __KRIG_sens__
#define __KRIG_sens__

#include <RcppEigen.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Combinatorial variance computation.
//' @description For a given combination this function computes the associated variance for the
//' variable enumerated by the combination values.
//' @param KF values of the kernel integral evaluations.
//' @param comb Combination.
//' @param X Points in the grid.
//' @param Gamma Cube with integral results.
//' @return Real value of sensitivity.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @seealso For a complete application you can check the documentation of \code{\link{Krigvar}}.
//' @references
//' \insertRef{Kanova:2013}{KRIG}
//' \insertRef{Aronszajn:ThRK}{KRIG}
//' @export
// [[Rcpp::export]]
double Krigidx( const Eigen::VectorXd& KF,
                const Eigen::VectorXd& comb,
                const Eigen::MatrixXd& X,
                const Eigen::ArrayXd& Gamma );

//--------------------------------------------------------------------------------------------------
//' @title Combinatorial variance computation.
//' @description Computation of variance
//' @param KF values of the kernel integral evaluations.
//' @param Gamma Cube with integral results.
//' @return Real value of sensitivity.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' library( KRIG ) 
//' options( stringsAsFactors = FALSE )
//'  
//' kernel_1<-function( x, y ) exp( -0.5*(x-y)^2)
//' kernel_2<-function( x, y ) exp( -0.7*(x-y)^2)
//' kernel_3<-function( x, y ) exp( -0.1*(x-y)^2)
//' 
//' Kernels<-data.frame( kernel = c( 'kernel_1', 'kernel_2', 'kernel_3' ), 
//'                      min = c( -1, -1, -2 ), 
//'                      max = c( 1, 1, 2 ),
//'                      n = c( 100, 100, 100 ) )
//'                      
//' n<-20
//' X<-matrix( c( seq( -1, 1, length.out = n ), 
//'               seq( -1, 1, length.out = n ),
//'               seq( -2, 2, length.out = n ) ), n, 3 )
//'               
//' KI<-list_integrate_kernel( Kernels, X )
//' GK<-Kanova( Kernels, KI, X )
//'     
//' f<-function( x ) x[1] + 30 * x[2] + 60 * x[3]
//' Func<-apply( X, 1, FUN = f )
//'     
//' KF<-solve( GK$Kanova + diag( 1e-8, n, n ), Func )
//'     
//' SbI<-NULL
//' for ( j in 1:3 ) {
//'   CB<-combn( 1:3, j )  
//'   for ( l in 1:ncol( CB ) ) {
//'     SbI<-c( SbI, Krigidx( KF, CB[,l], X, GK$Gamma ) )
//'     names(SbI)[length(SbI)]<-paste( 'C.', paste( CB[,l], collapse='.' ), sep = '' )
//'   }
//' }
//'   
//' Var<-Krigvar( KF, GK$Gamma )
//'     
//' SVar<-sum( SbI / Var )
//' 
//' @references
//' \insertRef{Kanova:2013}{KRIG}
//' \insertRef{Aronszajn:ThRK}{KRIG}
//' 
//' @export
// [[Rcpp::export]]
double Krigvar( const Eigen::VectorXd& KF, 
                const Eigen::ArrayXd& Gamma );

#endif
