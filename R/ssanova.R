#___________________________________________________________________________________________________
# Vector with integrated values
vect_integrate_kern<-function( x, k, a, b, method = 'Kronrod' ) {
  integra<-function( x, k, a, b ) integral( k, xmin = a, xmax = b, method = method, y = x )
  K<-sapply( x, FUN = integra, k, a, b )
  return( K )
}

#___________________________________________________________________________________________________
array_integrate_kern<-function( kernels, X ) {
  Y<-NULL
  for ( i in 1:ncol(X) ) { # i<-1
    Y<-cbind( Y, vect_integrate_kern( X[,i],  )
  }
  return( Y )
}

#___________________________________________________________________________________________________
kern_integral<-function( kernels ) {
  a<-NULL
  for ( i in 1:nrow(kernels) ) {
    a<-c( a, integral2( kernels[i,2], kernels[i,3], kernels[i,4], kernels[i,3], kernels[i,4] ) )
  }
  return( a )
}

