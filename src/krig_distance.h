#ifndef __KRIG_distance__
#define __KRIG_distance__

#include <RcppEigen.h>
#include <cmath>
#include <omp.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//' @title Generic weighted p-distance
//' @description Many isotropic models can be defined employing a distance. In particular this 
//' function implements a distance employing weights and different powers.
//' @param x First vector.
//' @param y Second vector.
//' @param w Weights for every coordinate in the vectors.
//' @param p Powers for every coordinate in the vectors. To define a true norm every coordinate
//' has to be greater than 1.
//' @return Real value of the weighted p-distance.
//' @author Pedro Guarderas \email{pedro.felipe.guarderas@@gmail.com}.
//' @examples
//' library(KRIG)
//' 
//' n<-1e5
//' x<-runif( n )
//' y<-runif( n )
//' p<-rep(2.5,n)
//' w<-runif(n)
//' 
//' kw<-weight_pow_dist( x, y, w, p )
//' nw<-sum( w * abs( x - y )^p )
//' 
//' nw == kw
//' @export
// [[Rcpp::export]]
double weight_pow_dist( const Eigen::VectorXd& x, 
                        const Eigen::VectorXd& y, 
                        const Eigen::VectorXd& w, 
                        const Eigen::VectorXd& p );
#endif
