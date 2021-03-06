% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{spherical_kernel}
\alias{spherical_kernel}
\title{Spherical kernel.}
\usage{
spherical_kernel(h, phi, theta)
}
\arguments{
\item{h}{distance variable.}

\item{phi}{amplitude parameter.}

\item{theta}{smoothness parameter.}
}
\value{
Real value.
}
\description{
Isotropic kernel commonly employed in geostatistic.
}
\examples{
h<-4.0
phi<-1.0
theta<-100.0
spherical_kernel( h, phi, theta )
}
\author{
Pedro Guarderas \email{pedro.felipe.guarderas@gmail.com}.
}
