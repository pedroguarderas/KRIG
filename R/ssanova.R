#___________________________________________________________________________________________________
# Vector with integrated values
vect_integrate_kern<-function( x, k, a, b, method = 'Kronrod' ) {
  integra<-function( x, k, a, b ) integral( k, xmin = a, xmax = b, method = method, y = x )
  K<-sapply( x, FUN = integra, k, a, b )
  return( K )
}

#___________________________________________________________________________________________________
gamma_integrate<-function( x, k, a, b, method = 'Kronrod' ) {
  K<-function( s, x, y ) return( k(x,s) * k(s,y) )
  integra<-function( x, k, a, b ) {
      return( integral( k, xmin = a, xmax = b, method = method, x = x[1], y = x[2] ) )
  }
  G<-apply( x, 1, FUN = integra, K, a, b )
  return( G )
}
  
#___________________________________________________________________________________________________
array_integrate_kern<-function( kernels, X ) {
  Y<-NULL
  for ( i in 1:ncol(X) ) { # i<-1
    Y<-cbind( Y, vect_integrate_kern( X[,i],  kernels[i,2], kernels[i,3], kernels[i,4] )  )
  }
  return( Y )
}

#___________________________________________________________________________________________________
kern_integral<-function( kernels ) {
  a<-NULL
  for ( i in 1:nrow(kernels) ) {
    a<-c( a, integral2( kernels[i,2], kernels[i,3], kernels[i,4], kernels[i,3], kernels[i,4] )$Q )
  }
  return( a )
}
