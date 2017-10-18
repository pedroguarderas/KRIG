# Example 2 ----------------------------------------------------------------------------------------
library(RKHSENS)
options( stringsAsFactors = FALSE )


kernel_1<-function( x, y ) exp( -0.5*(x-y)^2)
kernel_2<-function( x, y ) exp( -0.7*(x-y)^2)
kernel_3<-function( x, y ) exp( -0.1*(x-y)^2)

Kernels<-data.frame( kernel = c( 'kernel_1', 'kernel_2', 'kernel_3' ), 
                     min = c( 0, -1, -5 ), 
                     max = c( 1, 1, 5 ),
                     n = c( 500, 500, 500 ) )

n<-20
X<-matrix( c( seq( 0, 1, length.out = n ), 
              seq( -1, 1, length.out = n ),
              seq( -5, 5, length.out = n ) ), n, 3 )

KI<-RKHKernInteg( Kernels, X )

# Returns Gamma array and Kanova matrix
GK<-RKHAnova( Kernels, KI, X )

f<-function( x ) abs( x[1] + 30 * x[2] + 60 * x[3] )
F<-apply( X, 1, FUN = f )

Sobol<-SobolIndex( GK, X, F, 1:3 )
Sobol$S / Sobol$Var
