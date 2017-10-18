# Example 4 ----------------------------------------------------------------------------------------
# Gaussian process regression
library( RKHSENS )
library( rgl )

# Observed -----------------------------------------------------------------------------------------
m<-50
x<-seq( -3, 6, length.out = m )
y<-x
f<-function(x,y){
  return(exp(-x^2-y^2)+0.5*exp(-(x-2)^2-(y-2)^2) - 0.5*exp(-(x-4)^2-(y-0.25)^2)) 
}
z<-outer(x,y,f)

# Prediction ---------------------------------------------------------------------------------------
n<-30
x1<-runif( n, -2, 5 )
x2<-runif( n, -2, 5 )
X<-cbind(x1,x2)
Y<-as.matrix(expand.grid(x,y))

Z<-NULL
for(i in 1:n) Z<-c( Z, f(x1[i],x2[i]) )
Z<-as.matrix( Z, n, 1 )

# Kernel -------------------------------------------------------------------------------------------
s<-10.0
t<-1.1
Kern<-function( x, y ) return( RKHKerExp( sum( (x-y)^2 ), s, t ) )

# Gaussian process estimation ----------------------------------------------------------------------
krgs<-RKHEstimate( Z, X, Y, Kern )
W<-matrix( krgs$W, m, m )

# Plotting the results -----------------------------------------------------------------------------
persp3d( x, y, W, col='darkgreen', alpha = 1.0 )
persp3d( x, y, z, col='gold', alpha = 0.6, add = TRUE )
points3d( x1, x2, Z, size = 8.0, col = 'dodgerblue3', add = TRUE )
