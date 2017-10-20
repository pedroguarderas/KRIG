# Example 3 ----------------------------------------------------------------------------------------
# Gaussian process regression
library( RKHSENS )

# Observed -----------------------------------------------------------------------------------------
m<-100
x<-seq( -3, 6, length.out = m )
f<-function(x){
  return( 1 + exp(-x^2) + 0.5 * exp(-(x-2)^2 ) - 0.5 * exp( -(x-4)^2 ) )
}
z<-sapply( x, f )

# Prediction ---------------------------------------------------------------------------------------
n<-10
X<-matrix( runif( n, -3, 6 ), n, 1 )

m<-100
Y<-matrix( seq( -3, 6, length.out = m ), m, 1 )

Z<-matrix( sapply( X, f ), n, 1 ) + rnorm( n, 0, 0.02 )

# Kernel -------------------------------------------------------------------------------------------
s<-10.0
t<-1.1
w<-1.2
p<-2.1
Kern<-function( x, y ) {
  h<-RKHWeightPowDist( x, y, w, p )
  return( RKHKerExp( h, s, t ) )
}

# Gaussian process estimation ----------------------------------------------------------------------
K = RKHCov( X, X, Kern, TRUE );
k = RKHCov( Y, X, Kern );
S = diag( 10, n, n );
krgs<-RKHEstimate( Z, X, Y, K, k, S, 2, 1 )
W<-krgs$W

# Plotting the results -----------------------------------------------------------------------------
ymin<-min( z, W[,1], Z[,1] )
ymax<-max( z, W[,1], Z[,1] )
plot( x, z, type = 'l', lwd = 2, col='gold', ylim = c( ymin, ymax ) )
points( Y, W[,1], col='darkgreen', type = 'l', lwd = 2 )
points( X, Z[,1], col = 'dodgerblue3', pch = 16, cex = 1.2 )

# plot( W[,1] - z, pch = 16, cex = 1, col = 'dodgerblue3' )


