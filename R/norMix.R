####---- Normal Mixtures  "norMix" -------
####---- ~~~~~~~~~~~~~~~   ######  -------
#### Object-oriented  S/R - functions for dealing with 1D normal mixtures.
#### -------------------------------------------------------------------------
#### Author: Martin Maechler, 20 Mar 1997
#### -------------------------------------------------------------------------

norMix <- function(mu, sig2 = rep(1, m), w = NULL,
                   name = NULL, long.name = FALSE)
{
    ## Purpose: Constructor for 'norMix' (normal mixture) objects
    ## -------------------------------------------------------------------------
    ## Arguments: mu: vector of means;  sig2: vector of  variances  sigma^2
    ##		w : vector of weights (adding to 1) -- default: equal
    ##		name : name attribute; constructed from (mu,sig2) by default
    ##		long.name : logical used for default \code{name} construction
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Mar 97, 14:58

    if(!is.numeric(mu)) stop("'mu' must be numeric!")
    m <- length(mu)
    if(length(sig2) == 1) sig2 <- rep(sig2, m)
    if(length(sig2) != m || !is.numeric(sig2)|| any(sig2 <=0))
        stop("'sig2' = sigma^2  must be > 0 with same length as 'mu'")
    if(is.null(w))
        w <- rep(1/m, m)
    else {
        if(length(w) != m || !is.numeric(w) || any(w<0))
            stop("'w' must be >= 0  with same length as 'mu'")
        s <- sum(w)
        if(abs(s-1) > 10*.Machine$double.eps) w <- w / sum(w)
    }
    if(is.null(name)) {
        sformat <- function(v) sapply(v, format, digits=1)
        pPar <- function(pp) {
            if(long.name)
                paste("(",paste(sformat(pp),  collapse= ","),")", sep="")
            else
                paste(sformat(pp), collapse= "")
        }
        name <- paste("NM",format(m),".",
                      pPar(mu), "_", pPar(sig2), sep="")
    }
    structure(name = name, class = "norMix",
              .Data = cbind(mu = mu, sig2 = sig2, w = w))
}

is.norMix <- function(obj)
{
  ## Purpose: is 'obj' a "norMix", i.e.  Normal Mixture object ?
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:38
  inherits(obj, "norMix") &&
  (!is.null(w <- obj[,"w"])) &&
  is.numeric(w) && all(w >= 0) && abs(sum(w)-1) < 1000*.Machine$double.eps
}

m.norMix <- function(obj) nrow(obj) ##  Number of components of  normal mixture

mean.norMix <- function(x, ...)
{
  ## Purpose: Return "true mean", i.e., the expectation of  a normal mixture.
  if(!is.norMix(x)) stop("'x' must be a 'Normal Mixture' object!")
  x[,"w"] %*% x[,"mu"]
}

var.norMix <- function(x, ...)
{
  ## Purpose: 'true' Variance, i.e. E[(X- E[X])^2]  for X ~ normal mixture.
  if(!is.norMix(x)) stop("'x' must be a 'Normal Mixture' object!")
  w <- x[,"w"]
  mu <- w %*% x[,"mu"]
  w %*% (x[,"sig2"] + (x[,"mu"]-mu)^2)
}

print.norMix <- function(x, ...)
{
    ## Purpose: print method for  "norMix" objects (normal mixtures)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Mar 97, 10:02
    ox <- x
    has.nam <- !is.null(nam <- attr(x,"name"))
    cat("'Normal Mixture' object",
        if(has.nam) paste("\t ``", nam, "''", sep=''), "\n")
    if(has.nam) attr(x, "name") <- NULL
    cl <- class(x);  cl <- cl[ cl != "norMix"] #- the remaining classes
    class(x) <- if(length(cl)>0) cl ## else NULL
    NextMethod("print", ...)
    invisible(ox)
}

dnorMix <- function(obj, x = NULL, xlim = NULL, n = 511)
{
  ## Purpose: density evaluation for "norMix" objects (normal mixtures)
  ## -------------------------------------------------------------------------
  ## Arguments: obj: Normal Mixture object;  x:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.norMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  if(is.null(x)) {
    if(is.null(xlim)) { ##-- construct "reasonable" abscissa values
      xlim <- mean.norMix(obj) + c(-3,3)*sqrt(var.norMix(obj))
    }
    x <- seq(xlim[1], xlim[2], length = n)
  }
  m <- m.norMix(obj) #-- number of components
  y <- numeric(length(x))
  w <- obj[,"w"]; mu <- obj[,"mu"]; sd <- sqrt(obj[,"sig2"])
  for(j in 1:m)
    y <- y + w[j] * dnorm(x, mean = mu[j], sd = sd[j])
  list(x = x, y = y)
}


## This is based on rmultz2() from S-news by Alan Zaslavsky & Scott Chasalow
## see also /usr/local/app/R/R_local/src/combinat/R/rmultz2.R
## Arg.names like  rbinom();  returns  n x p matrix
                                        # ---> ~/R/MM/STATISTICS/multinom-Dist.R
rmultinom <- function(n, size, prob) {
    K <- length(prob) # #{classes}
    matrix(tabulate(sample(K, n * size, repl = TRUE, prob) + K * 0:(n - 1),
                    nbins = n * K),
           nrow = n, ncol = K, byrow = TRUE)
}

rnorMix <- function(n, obj)
{
    ## Purpose: Generate random numbers according to "norMix" object `obj'
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 27 Jun 2002, 16:03
    mu <- obj[,"mu"]
    sd <- sqrt(obj[,"sig2"])
    nj <- rmultinom(n=1, size = n, prob = obj[,"w"])
    ## Easy version: *round* the `nj' to the nearest integers
    ## this is *inaccurate* for small n !
    sample(unlist(sapply(seq(nj),
                         function(j) rnorm(nj[j], mean = mu[j], sd = sd[j]))))
}


plot.norMix <-
    function(x, type = "l", n = 511, xout = NULL, xlim = NULL,
             xlab = "x", ylab = "f(x)", main = attr(x,"name"), lwd = 1.4,
             p.norm = TRUE, p.h0 = TRUE,
             parNorm = list(col = 2, lty = 2, lwd = 0.4),
             parH0   = list(col = 3, lty = 3, lwd = 0.4),
             ...)
{
    ## Purpose: plot method for  "norMix" objects (normal mixtures)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Mar 1997
    d.o <- dnorMix(x, n=n, x = xout)
    if(p.norm)
        dn <- dnorm(d.o$x, mean = mean.norMix(x), sd = sqrt(var.norMix(x)))
    plot(d.o, type=type, xlim = xlim, ylim = c(0,max(d.o$y, if(p.norm) dn)),
         main = main, xlab = xlab, ylab = ylab, lwd = lwd, ...)
    if(p.norm)	do.call("lines",  c(list(x=d.o$x, y=dn), parNorm))
    if(p.h0)	do.call("abline", c(list(h = 0), parH0))
    invisible(x)
}

lines.norMix <-
    function(x, type = "l", n = 511, xout = NULL, lwd = 1.4,
             p.norm = FALSE, parNorm = list(col = 2, lty = 2, lwd = 0.4), ...)
{
    ## Purpose: lines method for "norMix" objects (normal mixtures)
    ## -------------------------------------------------------------
    ## Author: Martin Maechler, Date: 27 Jun 2002, 16:10
    xlim <- if(is.null(xout)) par("usr")[1:2] # else NULL
    d.o <- dnorMix(x, n = n, x = xout, xlim = xlim)
    lines(d.o, type=type, lwd = lwd, ...)
    if(p.norm) {
        dn <- dnorm(d.o$x, mean = mean.norMix(x), sd = sqrt(var.norMix(x)))
        do.call("lines", c(list(x=d.o$x, y=dn), parNorm))
    }
    invisible()
}


r.norMix <- function(obj, x = NULL, xlim = NULL, n = 511, xy.return = TRUE)
{
  ## Purpose: Compute r := f / f0; f = normal mixture; f0 = "best" normal approx
  ## Author : Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.norMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  d.o <- dnorMix(obj, x=x, xlim=xlim, n=n)
  dn  <- dnorm(d.o$x, mean = mean.norMix(obj), sd = sqrt(var.norMix(obj)))
  if(xy.return) list(x = d.o$x, y= d.o$y / dn, f0 = dn) else d.o$y / dn
}
