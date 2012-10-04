\name{dcomplnorm}
\alias{dcomplnorm}
\title{Composite lognormal pdf}
\description{Computes the pdf of the composite lognormal distribution}
\usage{
dcomplnorm(x, spec, sigma = 1, theta = 1, log=FALSE, ...)
}
\arguments{
  \item{x}{scale or vector of positive values at which the pdf needs to be computed}
  \item{sigma}{the value of sigma parameter of the lognormal distribution, must be positive}
  \item{theta}{the value of theta parameter, the cutoff point, must be positive}
  \item{spec}{the specific distribution with which the lognormal distribution should be composited with}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf values computed at \code{x}}
\references{S. Nadarajah, S. A. A. Bakar, CompLognormal: An R Package for Composite Lognormal Distributions, submitted}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
y=dcomplnorm(x,"exp",rate=1)}
\keyword{Composite lognormal}
