library( data.table )
library( RKHSENS )
library( rgl )
library(viridis)

load( 'RData/mining_data.RData' )

I<-seq( 1, nrow( Copper ), 1 )
n<-length(I)
X<-as.matrix( Copper[ I, list( x1, x2, x3 ) ] )

Z<-as.matrix( Copper[ I, list( Z ) ], n, 1 )

m<-c( 20, 20, 20 )
x1_lim<-c( min( X[,1] ), max( X[,1] ) )
x2_lim<-c( min( X[,2] ), max( X[,2] ) )
x3_lim<-c( min( X[,3] ), max( X[,3] ) )
Y1<-seq( x1_lim[1], x1_lim[2], length.out = m[1] )
Y2<-seq( x2_lim[1], x2_lim[2], length.out = m[2] )
Y3<-seq( x3_lim[1], x3_lim[2], length.out = m[3] )
Y<-expand.grid( Y1, Y2, Y3 )
Y<-as.matrix( Y )

s<-1
t<-100
Kern<-function( x, y ) {
  h<-sqrt( sum( ( x - y )^2 ) )
  return( RKHKerSqrExp( h, s, t ) )
}

K = RKHCov( X, X, Kern, TRUE ) + diag( 0.00001, n, n )
k = RKHCov( Y, X, Kern )

KRIG<-RKHEstimate( Z = Z,
                   K = K, 
                   k = k,
                   G = matrix( 0, 1, 1),
                   g = matrix( 0, 1, 1),
                   type = "ordinary", 
                   typeinv = 'syminv' )

nbcol<-200
color = rev(viridis(nbcol))
W<-matrix( KRIG$Z, m[1], m[2] )
zcol  = cut(W, nbcol)


# Design plot by layers
persp3d( Y1, Y2, W, color=color[zcol], alpha = 0.6 )
points3d( X[,1], X[,2], 0, color=color[zcol], size = 8.0, col = 'dodgerblue3', add = TRUE )

filled.contour( Y1, Y2, W, nlevels = 100, color = viridis )


gc()
