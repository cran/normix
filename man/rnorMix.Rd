\name{rnorMix}
\alias{rnorMix}
\title{Generate ``Normal Mixture'' Distributed Random Numbers}
\description{
  Generate \code{n} random numbers, distributed according to a normal
  mixture.
}
\usage{
rnorMix(n, obj)
}
\arguments{
  \item{n}{the number of random numbers desired.}
  \item{obj}{an object of class \code{norMix}.}
}
\details{
  For a mixture of \eqn{m}, i.e., \code{m.norMix(obj)}, components,
  generate the number in each component as multinomial, and then use
  \code{\link{rnorm}} for each.
}
\value{
  numeric vector of length \code{n}.
}
\seealso{\code{\link{dnorMix}} for the density, and \code{\link{norMix}}
  for the construction and further methods.
}
\examples{
x <- rnorMix(5000, MW.nm10)
hist(x)# you don't see the claw
plot(density(x), ylim = c(0,0.6),
     main = "Estim. and true 'MW.nm10' density")
lines(MW.nm10, col = "orange")
}
\keyword{distribution}
