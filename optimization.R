# conjugate_gradient<-function( A, b, x, n, e ) {
#   r<-A %*% x - b
#   d<- -r
#   alpha<- -as.numeric( ( t(r) %*% d ) / ( t(d) %*% A %*% d ) )
#   x<-x + alpha * d
#   k<-0
#   while( k < n && sqrt( t(r)%*%r ) > e ) {
#     r<- A %*% x - b
#     beta<-as.numeric( ( t(d) %*% A %*% r ) / ( t(d) %*% d ) )
#     d<- -r + beta * d
#     alpha<- -as.numeric( ( t(r) %*% d ) / ( t(d) %*% A %*% d ) )
#     x<-x + alpha * d
#     k<-k + 1
#   }
#   return( x )
# }

conjugate_gradient<-function( A, b, x, n, e ) {
  r<-b - A %*% x 
  p<-r
  g<-as.numeric( t(r) %*% r )
  k<-0
  while ( k < n & g > e ) {
    y<-A %*% p
    alpha<-g / as.numeric( t(y) %*% p )
    x<-x + alpha * p
    r<-r - alpha * y
    g0<-g
    g<-as.numeric( t(r) %*% r )
    beta<-g / g0
    p<-r + beta * p
    k<-k + 1
  }
  return( x )
}