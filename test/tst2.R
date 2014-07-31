#___________________________________________________________________________________________________
library(RKHSENS)
options( stringsAsFactors = FALSE )

#___________________________________________________________________________________________________
kernel_1<-function( x, y ) min(x,y)
kernel_2<-function( x, y ) exp(-x*y)

kernels<-data.frame( id = 1:2, kern = c( 'kernel_1', 'kernel_2' ), linf = c(0,-1), lsup = c(1,1) )
X<-data.frame( x1 = seq(0,1,length.out=10), x2 = seq(-1,1,length.out=10) )
I<-array_integrate_kern( kernels, X )
alpha<-kern_integral( kernels )
K<-eval_kernels( kernels, X, I, alpha )
