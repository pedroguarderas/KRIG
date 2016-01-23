library(rgl)
library(RKHSENS)

m<-40
x<-seq(-2,5,length.out = m)
y<-x
f<-function(x,y){
  return(exp(-x^2-y^2)+0.5*exp(-(x-2)^2-(y-2)^2) - 0.5*exp(-(x-4)^2-(y-0.25)^2)) 
}
z<-outer(x,y,f)

n<-60
x1<-runif( n, -2, 5 )
x2<-runif( n, -2, 5 )

X<-cbind(x1,x2)
x0<-as.matrix(expand.grid(x,y))
Z<-NULL
for(i in 1:n) Z<-c(Z,f(x1[i],x2[i]))

s<-10.0
t<-1.1

# kern<-function(x,y) ktriRKH( distRKH( x, y ), alpha = 5 )
kern<-function(x,y) kexpRKH( distRKH( x, y ), s, t )

Z<-as.matrix(Z)
X<-as.matrix(X)
x0<-as.matrix(x0)
Kriging<-krigingSimpleRKH( Z, X, x0, kern )
K<-Kriging$K
k0<-Kriging$k0


# L<-chol( K )
# J<-chol2inv( L )
# Z0<-k0 %*% J %*% Z
Z0<-k0 %*% solve( K, Z )

Z0<-matrix( Z0, m, m )


persp3d( x, y, Z0, col='darkgreen', alpha = 1.0 )
persp3d( x, y, z, col='gold', alpha = 0.6, add = TRUE )
points3d( x1, x2, Z, size = 8.0, col = 'purple', add = TRUE )
