\name{Ker_Exp}
\alias{Ker_Exp}

\title{
Ker_Exp
}
\description{
Exponential kernel
}
\usage{
Ker_Exp(x, y, s, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ First coordinate vector
}
  \item{y}{ Second coordinate vector with the same length of \code{x}
}
  \item{s}{ Applitude parameter for the exponential kernel
}
  \item{t}{ Smoothing parameter for the exponential kernel
}
}
\details{ The exponential kernel is defined like: \cr
\deqn{s^2 exp( -\frac{\sum_{i=1}^n (x_i - y_i)^2}{t^2} )}
}
\value{
Evaluation of the kernel for the given parameters
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Julio Andrade, Pedro Guarderas, Daniel Lagos, Andrés López, Paúl Recalde, Edison Salazar.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
\keyword{kernel}
