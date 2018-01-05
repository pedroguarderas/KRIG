# Load example data --------------------------------------------------------------------------------
library( data.table )
library( RKHSENS )
library( rgl )
library(viridis)

path<-'/home/aju/Documents/Lectures/spatial_statistics/data/'

data_file<-paste0( path, 'samples.dat' )

SG<-read.csv( file = data_file, header = FALSE, skip = 3, sep = ' ', dec = '.', 
              colClasses = c( 'numeric', 'numeric', 'numeric' ) )
SG<-as.data.table( SG )
SG[ , V1 := NULL ]
setnames( SG, c( 'x1', 'x2', 'Z' ) )                       

# Brenda mines
data_file<-paste0( path, 'BrendaMines.dat' )

BM<-read.csv( file = data_file, header = FALSE, skip = 4, sep = ',', dec = '.' )
BM<-as.data.table( BM )
BM[ , V1 := NULL ]
setnames( BM, c( 'a', 's', 'x1', 'x2', 'x3', 'Z1', 'Z2', 'C1', 'C2', 'C3' ) )   

# Setting kernel -----------------------------------------------------------------------------------
# n<-500
# I<-sample( 1:nrow( SG ), size = n, replace = FALSE )
I<-seq( 1, nrow( SG ), 20 )
n<-length(I)
X<-as.matrix( SG[ I, list( x1, x2 ) ] )

Z<-as.matrix( SG[ I, list( Z ) ], n, 1 )

m<-c( 80, 80 )
x1_lim<-c( min( X[,1] ), max( X[,1] ) )
x2_lim<-c( min( X[,2] ), max( X[,2] ) )
Y1<-seq( x1_lim[1], x1_lim[2], length.out = m[1] )
Y2<-seq( x2_lim[1], x2_lim[2], length.out = m[2] )
Y<-expand.grid( Y1, Y2 )
Y<-as.matrix( Y )

s<-1
t<-100
Kern<-function( x, y ) {
  h<-sqrt( sum( ( x - y )^2 ) )
  return( RKHKerSpher( h, s, t ) )
}

K = RKHCov( X, X, Kern, TRUE )
k = RKHCov( Y, X, Kern )

KRIG<-RKHEstimate( Z = Z,
                   K = K, 
                   k = k,
                   G = matrix( 0, 1, 1),
                   g = matrix( 0, 1, 1),
                   type = "simple", 
                   typeinv = 'syminv' )

nbcol<-200
color = rev(viridis(nbcol))
W<-matrix( KRIG$Z, m[1], m[2] )
zcol  = cut(W, nbcol)

persp3d( Y1, Y2, W, color=color[zcol], alpha = 0.6 )
points3d( X[,1], X[,2], 0, color=color[zcol], size = 8.0, col = 'dodgerblue3', add = TRUE )

filled.contour( Y1, Y2, W, nlevels = 100, color = viridis )


gc()
