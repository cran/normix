\name{r.norMix}
\alias{r.norMix}
\title{Ratio of Normal Mixture to Corresponding Normal}
\description{
  Compute \eqn{r(x) = f(x)/ f0(x)} where \eqn{f()} is a normal mixture
  density and \eqn{f0} the normal density with the same mean and
  variance as \eqn{f}.
}
\usage{
r.norMix(obj, x = NULL, xlim = NULL, n = 511, xy.return = TRUE)
}
\arguments{
  \item{obj}{an object of class \code{norMix}.}
  \item{x}{numeric vector with abscissa values where to evaluate the
    density. Default is constructed from \code{n} (and \code{xlim} if
    specified).}
  \item{xlim}{range of abscissa values, used if \code{x == NULL}.  By default,
    xlim taken as mean plus/minus 3 standard deviations of the normal mixture.}
  \item{n}{number of abscissa values to generate if \code{x} is not specified.}
  \item{xy.return}{logical indicating if the result should be a list or
    just a numeric vector, see below.}
}
\value{
  It depends on \code{xy.return}.  If it's false,  a numeric vector of
  the same length as \code{x}, if true (as per default), a list that can
  be plotted, with components
  \item{x}{abscissa values corresponding to argument \code{x}.}
  \item{y}{corresponding values \eqn{r(x)}.}
  \item{f0}{values of the moment matching normal density \eqn{f0(x)}.}
}
\note{The ratio function is used in certain semi-parametric density
  estimation methods (and theory).
}
\examples{
  d3 <- norMix(m = 5*(0:2), w = c(0.6, 0.3, 0.1))
  plot(d3)
  rd3 <- r.norMix(d3)
  str(rd3)
  stopifnot(rd3 $ y  == r.norMix(d3, xy.ret = FALSE))
  par(new = TRUE)
  plot(rd3, type = "l", col = 3, axes = FALSE, xlab = "", ylab="")
  axis(4, col.axis=3)
}
\keyword{distribution}
