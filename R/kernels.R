#___________________________________________________________________________________________________
# Algunos núcleos para kriging
Ker_Reg<-function(x,y,s) {
  return( sum( x * y ) + s^2 )
}

Ker_Tri<-function(x,y,c,a) {
  return( c * max( a - sqrt( sum( x - y ) ), 0 ) )
}

Ker_Exp<-function(x,y,s,t) {
  return( (s^2) * exp( -sum( (x - y)^2 ) / t^2 ) )
}

Ker_Gauss<-function(x,y,s,t) {
  return( exp( -sqrt( sum( (x - y)^2 ) ) / t ) * s^2 )
}

Ker_Matern<-function(x,y,v,s,t) {
  h<-sqrt( sum( ( x - y )^2 ) )
  return( (s^2) * ( ( ( 2 / t ) * sqrt(v) * h )^v  ) * 
            besselK( 2 * sqrt(v) * h / t, v ) / ( ( 2^(v-1) ) * gamma( v ) ) )
}