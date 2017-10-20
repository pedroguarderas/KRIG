# Example 5 ----------------------------------------------------------------------------------------
# Gaussian process regression
library( RKHSENS )
library( rgl )

# Observed -----------------------------------------------------------------------------------------
m<-50
x<-seq( -3, 6, length.out = m )
y<-x
f<-function(x,y){
  return(1 + exp(-x^2-y^2)+0.5*exp(-(x-2)^2-(y-2)^2) - 0.5*exp(-(x-4)^2-(y-0.25)^2)) 
}
z<-outer(x,y,f)

# Prediction ---------------------------------------------------------------------------------------
n<-40
x1<-runif( n, -3, 6 )
y1<-runif( n, -3, 6 )
X<-cbind(x1,y1)

m<-50
x0<-seq( -3, 6, length.out = m )
y0<-x0
Y<-as.matrix( expand.grid( x0, y0 ) )

Z<-NULL
for(i in 1:n) Z<-c( Z, f(x1[i],y1[i]) )
Z<-as.matrix( Z, n, 1 )

# Kernel -------------------------------------------------------------------------------------------
s<-10.0
t<-1.1
w<-c( 1.2, 1.5 )
p<-c( 2.3, 2.1 )
Kern<-function( x, y ) {
  h<-RKHWeightPowDist( x, y, w, p )
  return( RKHKerExp( h, s, t ) )
}

# Gaussian process estimation ----------------------------------------------------------------------
K = RKHCov( X, X, Kern, TRUE );
k = RKHCov( Y, X, Kern );
S = diag( 0, n, n );
krgs<-RKHEstimate( Z, X, Y, K, k, S, 1, 1 )
W<-matrix( krgs$W, m, m )

# Plotting the results -----------------------------------------------------------------------------
persp3d( x0, y0, W, col='darkgreen', alpha = 0.8 )
persp3d( x, y, z, col='gold', alpha = 0.6, add = TRUE )
points3d( x1, y1, Z, size = 8.0, col = 'dodgerblue3', add = TRUE )
