#___________________________________________________________________________________________________
# First test for the integration routines
library(RKHSENS)
options( stringsAsFactors = FALSE )

#___________________________________________________________________________________________________
kernel<-function( x, y ) min(x,y)
x<-seq(0,1,length.out=100)

kernels<-data.frame( id = 1:2, kern = c( 'kernel', 'kernel' ), linf = c(0,0), lsup = c(1,1) )
X<-data.frame( x1 = x, x2 = x )
Y<-array_integrate_kern( kernels, X )
alpha<-kern_integral( kernels )
