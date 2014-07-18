#___________________________________________________________________________________________________
library(pracma)

#___________________________________________________________________________________________________
kernel<-function( x, y ) min(x,y)
x<-seq(0,1,length.out=100)

#___________________________________________________________________________________________________
# Vector with integrated values
vect_integrate_kern<-function( x, k, a, b, method = 'Kronrod' ) {
  integra<-function( x, k, a, b ) integral( k, xmin = a, xmax = b, method = method, y = x )
  K<-sapply( x, FUN = integra, k, a, b )
  return( K )
}

#___________________________________________________________________________________________________
options( stringsAsFactors = FALSE )
kernels<-data.frame( id = 1:2, kern = c( 'kernel', 'kernel' ), linf = c(0,0), lsup = c(1,1) )
X<-data.frame( x1 = x, x2 = x )

array_integrate_kern<-function( kernels, X ) {
  Y<-NULL
  for ( i in 1:ncol(X) ) { # i<-1
    Y<-cbind( Y, vect_integrate_kern( X[,i], kernels[i,2], kernels[i,3], kernels[i,4] ) )
  }
  return( Y )
}

Y<-array_integrate_kern( kernels, X )

