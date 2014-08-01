#___________________________________________________________________________________________________
library(RKHSENS)
options( stringsAsFactors = FALSE )

#___________________________________________________________________________________________________
kernel_1<-function( x, y ) exp( -0.2*(x-y)^2)
kernel_2<-function( x, y ) exp( -0.2*(x-y)^2)

kernels<-data.frame( id = 1:2, kern = c( 'kernel_1', 'kernel_2' ), linf = c(0,-1), lsup = c(1,1) )
X<-data.frame( x1 = seq(0,1,length.out=10), x2 = seq(-1,1,length.out=10) )
I<-array_integrate_kern( kernels, X )
alpha<-kern_integral( kernels )

# Returns Gamma array and KANOVA matrix
GK<-eval_kernels( kernels, X, I, alpha )
C<-chol( GK$KANOVA )
C<-chol2inv( C )

# E<-eigen( GK$KANOVA )
# E$values
