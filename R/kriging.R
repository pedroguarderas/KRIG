#___________________________________________________________________________________________________
# Simple kriging
kriging_simple<-function( Z, X, x0, k ) {
  n<-nrow(X)
  m<-nrow(x0)
  k0<-Matrix( 0, n, m )
  K<-Matrix( 0, n, n )
  
  for ( i in 1:n ) { # i<-1
    for ( j in i:n ) {
      K[i,j]<-k( X[i,], X[j,] )
      if ( j > i ) {
        K[j,i]<-K[i,j]
      }
    }
    for ( l in 1:m ) {
      k0[i,l]<-k( x0[l,], X[i,] )
    }
  }
  L<-chol( K )
  J<-chol2inv( L )
  Z0<-t( k0 ) %*% J %*% Z
  return( list( Z0 = Z0, K = K, k0 = k0, L = L, J = J ) )
}

kriging_simple_cg<-function( Z, X, x0, k, l, n, e ) {
  K<-NULL
  k0<-NULL
  L<-NULL
  d<-dim( X )
  for ( i in 1:d[2] ) { # i<-1
    K<-rbind( K, apply( X, c(2), FUN = k, X[,i] ) ) 
    k0<-rbind( k0, apply( x0, c(2), FUN = k, X[,i] ) )
  }
  for ( i in 1:dim(x0)[2] ) {
    L<-cbind( L, conjugate_gradient( K, k0[,i], l, n, e ) )
  }
  Z0<-t( L ) %*% Z
  return( list( Z0 = Z0, K = K, k0 = k0, L = L ) )
}

#___________________________________________________________________________________________________
# Ordinary kriging
kriging_ordinary<-function( Z, X, x0, k ) {
  n<-nrow(X)
  m<-nrow(x0)
  k0<-Matrix( 0, n, m )
  K<-Matrix( 0, n, n )
  
  for ( i in 1:n ) { # i<-1
    for ( j in i:n ) {
      K[i,j]<-k( X[i,], X[j,] )
      if ( j > i ) {
        K[j,i]<-K[i,j]
      }
    }
    for ( l in 1:m ) {
      k0[i,l]<-k( x0[l,], X[i,] )
    }
  }

  L<-chol( K )
  J<-chol2inv( L )
  ones<-matrix( 1, d[2], 1 )
  u<-as.numeric( ( 1.0 / ( t( ones ) %*% J %*% ones ) ) * ( t( ones ) %*% J %*% Z ) )
  Z0<-u + t( k0 ) %*% J %*% ( Z - u )
  return( list( Z0 = Z0, K = K, k0 = k0, L = L, J = J ) )
}

kriging_ordinary_cg<-function( Z, U, X, x0, u0, k, l, n, e ) {
  K<-NULL
  k0<-NULL
  L<-NULL
  d<-dim( X )
  for ( i in 1:d[2] ) { # i<-1
    K<-rbind( K, apply( X, c(2), FUN = k, X[,i] ) ) 
    k0<-rbind( k0, apply( x0, c(2), FUN = k, X[,i] ) )
  }
  
  for ( i in 1:dim(x0)[2] ) {
    L<-cbind( L, conjugate_gradient( K, k0[,i], l, n, e ) )
  }
  Z0<-u0 + t( L ) %*% ( Z - U )
  return( list( Z0 = Z0, K = K, k0 = k0, L = L ) )
}
