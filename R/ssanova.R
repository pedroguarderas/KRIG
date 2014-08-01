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

#___________________________________________________________________________________________________
# Evaluate kernel
eval_kernels<-function( kernels, X, I, alpha ) {
  n<-nrow(X)
  m<-ncol(X)
  Gamma<-array( 0, dim = c( n, n, m ) )
  KANOVA<-Matrix( 1, n, n )
  for ( i in 1:m ) {
    for ( k in 1:n ) {
      for ( l in k:n ) {
        Gamma[k,l,i]<-eval( call( kernels[i,2],  X[k,i], X[l,i] ) ) - I[k,i] * I[l,i] / alpha[i] 
        KANOVA[k,l]<-KANOVA[k,l] * ( 1 + Gamma[k,l,i] )
        if ( l > k ) {
          Gamma[l,k,i]<-Gamma[k,l,i]
          KANOVA[l,k]<-KANOVA[k,l]
        }
      }
    }
  }
  return( list( Gamma = Gamma, KANOVA = KANOVA ) )
}

#___________________________________________________________________________________________________
# Sobol indices
SobolIndex<-function( Eval, X, F ) {
  n<-nrow(X)
  m<-ncol(X)
  
  U<-solve( Eval$KANOVA, F )

  V<-Matrix( 1, n, n )
  for( i in 1:m ) {
    V<-V * ( Matrix( 1, n, n ) + Eval$Gamma[,,i] )
  }
  V<-V - Matrix( 1, n, n )
  
  S<-NULL
  for ( j in 1:m ) { # j<-1
    CMB<-combn( 1:m, j )
    for ( k in 1:ncol(CMB) ) {
      W<-Matrix( 1, n, n )
      for ( l in 1:nrow(CMB) ) {
        W<-W * Eval$Gamma[,,CMB[l,k]]
      }
      S<-c( S, as.numeric( t(U) %*% W %*% U ) )
      names(S)[length(S)]<-paste( CMB[,k], collapse='' )
    }
  }
  Var<-as.numeric( t(U) %*% V %*% U ) 
  return( list( S = S, Var = Var ) )
}




