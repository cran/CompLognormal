\name{qcomplnorm}
\alias{qcomplnorm}
\title{Composite lognormal quantile}
\description{Computes the quantile function of the composite lognormal distribution}
\usage{
qcomplnorm(p, spec, sigma = 1, theta = 1, log.p=FALSE, ...)
}
\arguments{
  \item{p}{scale or vector of probabilities at which the quantile function needs to be computed}
  \item{sigma}{the value of sigma parameter of the lognormal distribution, must be positive}
  \item{theta}{the value of theta parameter, the cutoff point, must be positive}
  \item{spec}{the specific distribution with which the lognormal distribution should be composited with}
  \item{log.p}{if TRUE then quantiles are returned for exp(p)}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{p}, giving the quantile values computed at \code{p}}
\references{S. Nadarajah, S. A. A. Bakar, CompLognormal: An R Package for Composite Lognormal Distributions, submitted}
\author{Saralees Nadarajah}
\examples{p=runif(10,min=0,max=1)
y=qcomplnorm(p,"exp",rate=1)}
\keyword{Composite lognormal}
