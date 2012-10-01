\name{nlogl.complnorm}
\alias{nlogl.complnorm}
\title{Composite lognormal likelihood}
\description{Computes the negative log-likelihood for a random sample from the composite lognormal distribution}
\usage{
nlogl.complnorm(p, x, spec, ...)
}
\arguments{
  \item{x}{the observed data}
  \item{p}{a vector of parameters of the composite lognormal distribution (the first and second parameters should be sigma and theta respectively}
  \item{spec}{the specific distribution with which lognormal should be composited with}
  \item{...}{other parameters}
}
\value{Negative log-likelihood function at the parameter values specified by \code{p} for the composite lognormal distribution}
\references{S. Nadarajah, S. A. A. Bakar, CompLognormal: A new R package, submitted}
\author{Saralees Nadarajah}
\examples{x=rexp(100)
nlogl.complnorm(p=c(1,1,1),x,"exp",rate)}
\keyword{Composite lognormal}
