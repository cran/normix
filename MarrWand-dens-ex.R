## Load Martin Maechler's   'Normal Mixtures' package :
source("/u/maechler/R/MM/STATISTICS/NorMix.S")

## Load the Marron & Wand examples:
source("/u/maechler/R/MM/STATISTICS/MarrWand-dens.S")

##--- Test the  'NorMix' package with these:

sapply(1:16, function(i) is.NorMix(get(paste("MW.nm",i,sep = ""))))
## Nice summary :
for(n in 1:16) {
    nam <- paste("MW.nm",n,sep = "")
    cat(nam,":\n")
    print(get(nam))
}

round(sapply(1:16, function(i) mean.NorMix(get(paste("MW.nm",i,sep = "")))),7)
      sapply(1:16, function(i) var.NorMix(get(paste("MW.nm",i,sep = ""))))

mult.fig(16, main = "Marron Wand densities")
for(n in 1:16)
  plot(get(paste("MW.nm",n,sep = "")), n = 1 + if(any(n == 14:15)) 1024 else 4096)

mult.fig(16, main = expression(r(x) == f(x) / f[0](x) *
    "  for Marron Wand densities"), line.main = 0.5)
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  plot(r.NorMix(nm), xlab = "x", ylab = "r(x)", type = 'l', main = attr(nm,"name"))
}

if(is.R()) require(modreg) # smooth.spline

mult.fig(16, main = "Marron Wand -- r'' ;  r = f / f0 --- Sm.Spl df=80")
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  sp <- smooth.spline(r <- r.NorMix(nm, n = 1001), df = 80)
  plot(predict(sp, deriv = 2),
       xlab = "x", ylab = "r''(x)", type = 'l', main = attr(nm,"name"))
}

mult.fig(16, main = "Marron Wand -- r'' ;  r = f / f0 --- Sm.Spl df=90")
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  sp <- smooth.spline(r <- r.NorMix(nm, n = 8001), df = 90)
  plot(predict(sp, deriv = 2),
       xlab = "x", ylab = "r''(x)", type = 'l', main = attr(nm,"name"))
}

mult.fig(16, main = "Marron Wand -- r'' ;  r = f / f0 --- Sm.Spl   GCV")
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  sp <- smooth.spline(r <- r.NorMix(nm, n = 4001))
  plot(predict(sp, deriv = 2), sub = paste("equiv.df=",format(sp$df)),
       xlab = "x", ylab = "r''(x)", type = 'l', main = attr(nm,"name"))
}


d2.est <- function(y, r, dx = 1)
{
  ## Purpose: Compute "optimal" estimate for f''(x_i), given y_i = f(x_i) + Z_i
  ## ~~~~~~~  assuming EQUIdistant x_i,  i=1,..,n
  ## -------------------------------------------------------------------------
  ## Arguments: y: vector (y_i), assuming  y_i = f(x_i) + Z_i,  i=1,..n
  ##            r: "half-bandwidth": integer >=1; SMOOTHING param (see below)
  ##	       dx: == x_{i+1} - x_i  (must be the same for all i !)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 1990 - Jan 93
  ## Ref.:   Ph.D. thesis of  Martin B. M"achler, Sem.f.Stat., ETH Z"urich
  ## -------------------------------------------------------------------------
  r <- as.integer(r)
  if(r < 1)  stop("'r' (half-bandwidth) must be integer >= 1")
  ##convolveFFT(y, wk21(r)/(dx^2))
  convolve(y, wk21(r)/(dx^2), type = "filter")
}

wk21 <- function(r)
{
  ## Purpose: Weights for  Smoothing f''(x)  for EQUIDISTANT x --> my thesis
  ## -------------------------------------------------------------------------
  ## Arguments: r: integer = half-bandwidth, rr1 = r(r+1)  is (odd) bandwidth.
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 1988
  rr1 <- r * (1 + r)
  b <- 1 + 2 * r
  k <- seq( - r, r)
  30 * (3 * k^2 - rr1)/(rr1 * b * (b-2) * (b+2))
}

np <- 4001; r <- 50
mult.fig(16, main = paste("Marron Wand -- r'' ;  r = f / f0 --- N=",
		 np,";  d2(r=", r,")"))
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  rN <- r.NorMix(nm, n = np);  x <- rN$x
  plot(x[-c(1:r,np+1-(r:1)) ], d2.est(rN$y,dx = x[2]-x[1], r = r),
       xlab = "x", ylab = "r''(x)", type = 'l', main = attr(nm,"name"))
}

np <- 4001; r <- 10 ; iix <- -c(1:r,np+1-(r:1))
mult.fig(16, main = paste("Marron Wand -- r'' * f0 ;  r = f / f0 --- N=",
		 np,";  d2(r=", r,")"))
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  rN <- r.NorMix(nm, n = np);  x <- rN$x
  plot(x[iix], rN$f0[iix] * d2.est(rN$y,dx = x[2]-x[1], r = r),
       xlab = "x", ylab = "r''(x) * f0", type = 'l', main = attr(nm,"name"))
}


if(is.R()) require(SfS) # for D1() :

###-------- MM's YASPE is bad..... :
np <- 4001
mult.fig(16, main = paste("Marron Wand -- r = f / f0  --- (r'/f0)' --- N=",np))
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  rN <- r.NorMix(nm, n = np);  x <- rN$x
  d1r <- D1(x, rN$y)
  r1f0 <- d1r / rN$f0
  plot(predict(smooth.spline(x,r1f0), deriv = 1),
       xlab = "x", ylab = "(r'/f0)'(x)", type = 'l', main = attr(nm,"name"))
}


np <- 1001
mult.fig(16, main = paste("Marron Wand -- b(x) ^2  Hjort-Glad  vs. YASPE"))
for(n in 1:16) {
  nm <- get(paste("MW.nm",n,sep = ""))
  rN <- r.NorMix(nm, n = np);  x <- rN$x
  d1r <- D1(x, rN$y)
  r1f0 <- d1r / rN$f0
  ## (r' / f0)' :
  b.yaspe <- predict(smooth.spline(x,r1f0), deriv = 1)$y
  ## f0*(r')' = f0 r''
  b.HjoGl <- rN$f0 * predict(smooth.spline(x,d1r), deriv = 1)$y
  yl <- pmax(1e-12,range((b.y2 <- b.yaspe^2),(b.HG2 <- b.HjoGl^2)))
  plot(x, b.y2, ylim = yl, xlab = "x", ylab = "b^2 (x)",
       type = "l", log = "y", err = -1)
  lines(x, b.HG2, col = 2, err = -1)
  par(new = TRUE)
  plot(nm, xx = x, p.norm = FALSE, p.h0 = FALSE, axes = FALSE,
       col = 3, lty = 2, lwd = 0.1)
}
