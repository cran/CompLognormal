\name{rcomplnorm}
\alias{rcomplnorm}
\title{Composite lognormal random numbers}
\description{Generates random numbers from the composite lognormal distribution}
\usage{
rcomplnorm(n, spec, sigma = 1, theta = 1, ...)
}
\arguments{
  \item{n}{number of random numbers to be generated}
  \item{sigma}{the value of sigma parameter of the lognormal distribution, must be positive}
  \item{theta}{the value of theta parameter of the lognormal distribution, must be positive}
  \item{spec}{the specific distribution with which lognormal should be composited with}
  \item{...}{other parameters}
}
\value{An object of the length \code{n}, giving the random numbers from the composite lognormal distribution}
\references{S. Nadarajah, S. A. A. Bakar, CompLognormal: A new R package, submitted}
\author{Saralees Nadarajah}
\examples{y=rcomplnorm(100,"exp",rate=1)}
\keyword{Composite lognormal}
