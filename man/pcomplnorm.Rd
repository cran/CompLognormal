\name{pcomplnorm}
\alias{pcomplnorm}
\title{Composite lognormal cdf}
\description{Computes the cdf of the composite lognormal distribution}
\usage{
pcomplnorm(x, spec, sigma = 1, theta = 1, log.p=FALSE, lower.tail=TRUE, ...)
}
\arguments{
  \item{x}{scale or vector of positive values at which the cdf needs to be computed}
  \item{sigma}{the value of sigma parameter of the lognormal distribution, must be positive}
  \item{theta}{the value of theta parameter, the cutoff point, must be positive}
  \item{spec}{the specific distribution with which the lognormal distribution should be composited with}
  \item{log.p}{if TRUE then log(cdf) are returned}
  \item{lower.tail}{if TRUE then cdf are returned else 1-cdf are retured}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the cdf values computed at \code{x}}
\references{S. Nadarajah, S. A. A. Bakar, CompLognormal: An R Package for Composite Lognormal Distributions, submitted}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
y=pcomplnorm(x,"exp",rate=1)}
\keyword{Composite lognormal}
