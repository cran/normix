\name{rmultinom}
\alias{rmultinom}
\title{Multinomial Distribution  Random Numbers}
%%---- FIXME: add this to base R -- together with "dmultinom"
\description{
  Generate multinomially distributed random number vectors.
}
\usage{
rmultinom(n, size, prob)
}
\arguments{
  \item{n}{number of random vectors to draw.}
  \item{size}{integer, say \eqn{N}, specifying the total number
   of objects that are put into \eqn{K} boxes in the typical multinomial
   experiment.}%% fixme
  \item{prob}{numeric non-negative vector of length \eqn{K}, specifying
    the probability for the \eqn{K} classes; is internally normalized to
    sum 1.}
}
\details{
  By definition, each component \code{j} is binomially distributed as
  \code{Bin(size, prob[j])}.

  The current implementation is based on \code{\link{sample}} and
  \code{\link{tabulate}} and hence a bit slower than optimally.
}
\value{
  an integer \code{n x K} matrix where each column is a random vector
  generated according to the desired multinomial law, and hence summing
  to \code{size}.
}
\references{The original idea was posted on the S-news mailing
  by Alan Zaslavsky and later modified by Scott Chasalow.
}
\seealso{\code{\link{rbinom}} which is a special case conceptually.
%% but does not return 2-vectors
}

\examples{
rmultinom(10, size = 12, prob=c(0.1,0.2,0.8))

N <- 20
pr <- c(1,3,6,10)
pr <- pr/sum(pr) # << normalization not necessary for generation

set.seed(153)
rr <- rmultinom(5000, N, prob = pr)
stopifnot(apply(rr,1,sum) == N)
(m <- apply(rr, 2, mean))
all.equal(m, N * pr/sum(pr)) # rel.error 0.00324
}
\keyword{distribution}
