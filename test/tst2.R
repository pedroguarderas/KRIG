#___________________________________________________________________________________________________
library(RKHSENS)
options( stringsAsFactors = FALSE )

#___________________________________________________________________________________________________
kernel_1<-function( x, y ) exp( -0.5*(x-y)^2)
kernel_2<-function( x, y ) exp( -0.7*(x-y)^2)
kernel_3<-function( x, y ) exp( -0.1*(x-y)^2)

kernels<-data.frame( id = 1:3, 
                     kern = c( 'kernel_1', 'kernel_2', 'kernel_3' ), 
                     linf = c(0,-1, 0 ), 
                     lsup = c(1,1,2) )
X<-data.frame( x1 = runif( 30, 0, 1 ), 
               x2 = runif( 30, -1, 1 ),
               x3 = runif( 30, 0, 2 ))
I<-array_integrate_kern( kernels, X )
alpha<-kern_integral( kernels )

# Returns Gamma array and KANOVA matrix
GK<-eval_kernels( kernels, X, I, alpha )

f<-function( x ) abs( x[1] + 30 * x[2] + 60 * x[3] )
F<-apply( X, 1, FUN = f )

Sobol<-SobolIndex( GK, X, F )
Sobol$S / Sobol$Var
