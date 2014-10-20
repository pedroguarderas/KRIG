\name{kriging_ordinary}
\alias{kriging_ordinary}

\title{
kriging_ordinary
}
\description{
Function for data interpolation by kriging, it uses a \code{ordinary Kriging}. Also computes
the mean.
}
\usage{
kriging_ordinary(Z, X, x0, k)
}
\arguments{
  \item{Z}{
  Data matrix for estimation point.
}
  \item{X}{
  Data matrix of known values of the function in the neighborhood of the point.
}
  \item{x0}{
  Matrix of correlations for the point.
}
  \item{k}{
  Function \code{kernel} to be applied in the interpolation by kriging.
}
}
\details{
Function \code{kriging_ordinary} is useful in the case of prediction with the known mean.
}
\value{
Function \code{kriging_ordinary} returns an object of class \code{list} containing the following components:
\item{Z0}{Matrix of observed values.}
\item{K}{Matrix of correlation.}
\item{L}{Covariance Matrix positive-definite square from Choleski factorization.}
\item{J}{Inverted Covariance Matrix \code{L}.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Julio Andrade, Pedro Guarderas, Daniel Lagos, Andrés López, Paúl Recalde, Edison Salazar.
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{kriging_ordinary_cg}}, \code{\link{Ker_Exp}}
}
\examples{
library(RKHSENS)
library(rgl)

m<-40
x<-seq(-2,5,length.out=m)
y<-x
f<-function(x,y){
  return(1+exp(-x^2-y^2)+0.5*exp(-(x-2)^2-(y-2)^2) - 0.5*exp(-(x-4)^2-(y-0.25)^2)) 
}
z<-outer(x,y,f)

n<-30
x1<-runif( n, -2, 5 )
x2<-runif( n, -2, 5 )

X<-cbind(x1,x2)
x0<-as.matrix(expand.grid(x,y))
Z<-NULL
for(i in 1:n) Z<-c(Z,f(x1[i],x2[i]))

s<-10.0
t<-1.1
k<-function( x, y ) return( Ker_Exp( x, y, s, t ) )
krgs<-kriging_ordinary( Z, X, x0, k )
Z0<-matrix( krgs$Z0, m, m )


persp3d( x, y, Z0, col='darkgreen', alpha = 1.0 )
persp3d( x, y, z, col='gold', alpha = 0.6, add = TRUE )
points3d( x1, x2, Z, size = 8.0, col = 'purple', add = TRUE )
}
\keyword{kriging}
