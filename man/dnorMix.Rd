\name{dnorMix}
\alias{dnorMix}
\title{Normal Mixture Density}
\description{
  Evaluate the density function of the normal mixture specified as
  \code{\link{norMix}} object.
}
\usage{
dnorMix(obj, x = NULL, xlim = NULL, n = 511)
}
\arguments{
  \item{obj}{an object of class \code{norMix}.}
  \item{x}{numeric vector with abscissa values where to evaluate the
    density. Default is constructed from \code{n} (and \code{xlim} if
    specified).}
  \item{xlim}{range of abscissa values, used if \code{x == NULL}.  By default,
    xlim taken as mean plus/minus 3 standard deviations of the normal mixture.}
  \item{n}{number of abscissa values to generate if \code{x} is not specified.}
}
\value{
  A list with components
  \item{x}{the abscissa values.}
  \item{y}{the density \eqn{f(x)} at \code{x}.}
}
\seealso{\code{\link{rnorMix}} for random number generation, and
  \code{\link{norMix}} for the construction and further methods,
  particularly \code{\link{plot.norMix}} which makes use \code{dnorMix}.
}
\examples{
 ff <- dnorMix(MW.nm7)
 str(ff)
 plot(ff, type = "h") # rather use plot(ff, ...)
}
\keyword{distribution}
