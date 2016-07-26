#___________________________________________________________________________________________________
library(RKHSENS)
options( stringsAsFactors = FALSE )

#___________________________________________________________________________________________________
N<-10
kernel_1<-function( x, y ) exp( -0.5*(x-y)^2)
kernel_2<-function( x, y ) exp( -0.7*(x-y)^2)
kernel_3<-function( x, y ) exp( -0.1*(x-y)^2)

kernels<-data.frame( id = 1:3, 
                     kern = c( 'kernel_1', 'kernel_2', 'kernel_3' ), 
                     linf = c( 0, -1, -5 ), 
                     lsup = c( 1, 1, 5 ) )
X<-data.frame( x1 = seq( 0, 1, length.out = N ), 
               x2 = seq( -1, 1, length.out = N ),
               x3 = seq( -5, 5, length.out = N ) )
I<-array_integrate_kern( kernels, X )
alpha<-kern_integral( kernels )

# Returns Gamma array and KANOVA matrixB
system.time( GK<-eval_kernels( kernels, X, I, alpha ) )

f<-function( x ) abs( x[1] + 30 * x[2] + 60 * x[3] )
F<-apply( X, 1, FUN = f )

Sobol<-SobolIndex( GK, X, F, 1:3 )
Sobol$S / Sobol$Var

#___________________________________________________________________________________________________
Kernels<-list()
Kernels[[1]]<-list( kernel = kernel_1, grid = seq( 0, 1, length.out = N ), min =  0, max = 1, size = 100 )
Kernels[[2]]<-list( kernel = kernel_2, grid = seq( -1, 1, length.out = N ), min = -1, max = 1, size = 100 )
Kernels[[3]]<-list( kernel = kernel_3, grid = seq( -5, 5, length.out = N ), min = -5, max = 5, size = 100 )

Integral<-integralRKH( Kernels )
system.time( KerEval<-evalKernRKH( Kernels, Integral ) )


