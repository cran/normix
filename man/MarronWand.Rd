\name{MarronWand}
\alias{MarronWand}% to use as "link"
%% Note: These are ``datasets'' , but cannot be handled as strictly
%%       as such since they  *require* the normix package.
\alias{MW.nm1}
\alias{MW.nm2}
\alias{MW.nm3}
\alias{MW.nm4}
\alias{MW.nm5}
\alias{MW.nm6}
\alias{MW.nm7}
\alias{MW.nm8}
\alias{MW.nm9}
\alias{MW.nm10}
\alias{MW.nm11}
\alias{MW.nm12}
\alias{MW.nm13}
\alias{MW.nm14}
\alias{MW.nm15}
\alias{MW.nm16}
\non_function{}
\title{Marron-Wand Densities as `norMix' Objects}
\description{
  The fifteen density examples used in Marron and Wand (1992)'s
  simulation study have been used in quite a few subsequent studies,
  can all be written as normal mixtures and provided here for
  convenience and didactical examples of normal mixtures.
  Number 16 has been added by Jansen et al.
}
\usage{
 MW.nm1
 MW.nm2
 MW.nm3
 MW.nm4
 MW.nm5
 MW.nm6
 MW.nm7
 MW.nm8
 MW.nm9
 MW.nm10
 MW.nm11
 MW.nm12
 MW.nm13
 MW.nm14
 MW.nm15
 MW.nm16
}
\source{
  The first part % 1..15 ??
  are translated from Steve Marron's Matlab code at
  \url{ftp://ftp.stat.unc.edu/pub/papers/marron/parameters/nmpar.m}
}
\author{Martin Maechler}
\references{
  Marron, S. and Wand, M. (1992)
  .............. %% FIXME

  For number 16,
  Janssen, et. al. (...) %% FIXME
}
\examples{
MW.nm10
plot(MW.nm14)

## How are they defined:
## 1.  as norMix() calls in  ../R/zMarrWand-dens.R
## 2.
require(normix)
ppos <- which("package:normix" == search())
nms <- ls(pat="^MW.nm", pos = ppos)
nms <- nms[order(as.numeric(substring(nms,6)))]
for(n in nms) {
   cat("\n",n,":\n"); print(get(n, pos = ppos))
}

## Plot all of them:
op <- par(mfrow=c(4,4), mgp = c(1.2, 0.5, 0), tcl = -0.2,
          mar = .1 + c(2,2,2,1), oma = c(0,0,3,0))
for(n in nms) { plot(get(n, pos = ppos))}
mtext("The Marron-Wand Densities", outer = TRUE,
      font = 2, cex = 1.6)
par(op)
}
\keyword{datasets}
\keyword{distribution}
