\name{norMix}
\alias{norMix}
\alias{is.norMix}
\alias{m.norMix}
\alias{mean.norMix}
\alias{var.norMix}
\alias{print.norMix}
\title{Mixtures of Univariate Normal Distributions}
\description{
  Objects of class \code{norMix} represent finite mixtures of
  (univariate) normal (aka Gaussian) distributions.  Methods for
  construction, printing, plotting, and basic computations are provided.
}
\usage{
norMix(mu, sig2 = NULL, w = NULL, name = NULL, long.name = FALSE)

is.norMix(obj)
m.norMix(obj)
mean.norMix(x, \dots)
var.norMix(x, \dots)
\method{print}{norMix}(x, \dots)
}
\arguments{
  \item{mu}{numeric vector of length \eqn{K}, say, specifying the means
    \eqn{\mu} of the \eqn{K} normal components.}
  \item{sig2}{numeric vector of length \eqn{K}, specifying the variances
    \eqn{\sigma^2} of the \eqn{K} normal components.}
  \item{w}{numeric vector of length \eqn{K}, specifying the mixture proportions
    \eqn{\pi_j}{p[j]} of the normal components, \eqn{j = 1,\dots,K}.
    Defaults to equal proportions}
  \item{name}{optional name tag of the result (used for printing).}
  \item{long.name}{ ~~Describe \code{long.name} here~~ }

  \item{obj,x }{an object of class \code{norMix}.}
  \item{\dots}{further arguments passed to methods.}
}
\details{
  The (one dimensional) normal mixtures, \R objects of class
  \code{"norMix"}, are constructed by \code{norMix} and tested for by
  \code{is.norMix}.  \code{m.norMix()} returns the number of mixture
  components,
  ....
 %% FIXME
  For further methods see below. %% seealso
}
\value{
  \code{norMix} returns objects of class \code{"norMix"} which are
  currently implemented as 3-column matrix with column names \code{mu},
  \code{sig2}, and \code{w}, and further attributes.
  The user should rarely need to access the underlying structure
  directly.
}
\author{Martin Maechler}
\seealso{\code{\link{dnorMix}} for the density, \code{\link{rnorMix}}
  for random numbers and \code{\link{print.norMix}}, the plot method.

  \code{\link{MarronWand}} has the Marron-Wand densities as normal mixtures.
}
\examples{
ex <- norMix(mu = c(1,2,5))# s^2 = 1, equal proportions
ex
plot(ex)# looks like a mixture of only 2
}
\keyword{distribution}
