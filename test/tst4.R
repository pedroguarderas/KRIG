#___________________________________________________________________________________________________
# Kriging example in 3D
library(rgl)
library(RKHSENS)

m<-40
x<-seq(-2,5,length.out=m)
y<-x
f<-function(x,y){
  return(exp(-x^2-y^2)+0.5*exp(-(x-2)^2-(y-2)^2) - 0.5*exp(-(x-4)^2-(y-0.25)^2)) 
}
z<-outer(x,y,f)

n<-50
x1<-runif( n, -2, 5 )
x2<-runif( n, -2, 5 )

X<-cbind(x1,x2)
x0<-as.matrix(expand.grid(x,y))
Z<-NULL
for(i in 1:n) Z<-c(Z,f(x1[i],x2[i]))

s<-10.0
t<-1.1
k<-function( x, y ) return( Ker_Exp( x, y, s, t ) )
krgs<-kriging_simple( Z, X, x0, k )
Z0<-matrix( krgs$Z0, m, m )


persp3d( x, y, Z0, col='dodgerblue', alpha = 1.0 )
persp3d( x, y, z, col='gold', alpha = 0.6, add = TRUE )
points3d( x1, x2, Z, size = 8.0, col = 'purple', add = TRUE )

