#### {filename starting with 'z' : must come after other the Normix() defining code}

### The first part of this is from  Steve Marron,
### from his Matlab code
###         ftp://ftp.stat.unc.edu/pub/papers/marron/parameters/nmpar.m
##
##    For generating Normal mixture parameters,
##      of the 15 Normal Mixture densities
##      from Marron and Wand (1992),
##      plus number 16 from Janssen, et. al.

MW.nm1 <- norMix(name = "#1 Gaussian",
		 mu = 0, sig2 = 1, w = 1)

MW.nm2 <- norMix(name = "#2 Skewed",
		 mu   = c(-.3, .3, 1),
		 sig2 = c(1.44, .64, 4/9),
		 w    = c(.2, .2, .6))

sig <- (2/3)^(0:7)
MW.nm3 <- norMix(name = "#3 Str Skew",
		 mu   = 3 * (sig - 1),
		 sig2 = sig^2,
		 w    = rep(1,8) / 8)

MW.nm4 <- norMix(name = "#4 Kurtotic",
		 mu   = c(0, 0),
		 sig2 = c(1 , .01),
		 w    = c(2/3, 1/3))


MW.nm5 <- norMix(name = "#5 Outlier",
		 mu   = c(0, 0),
		 sig2 = c(1, .01),
		 w    = c(.1, .9))

MW.nm6 <- norMix(name = "#6 Bimodal",
		 mu   = c(-1, 1),
		 sig2 = c(4/9, 4/9),
		 w = c(.5, .5))

MW.nm7 <- norMix(name = "#7 Separated",
		 mu  = c(-1.5, 1.5),
		 sig2 = c(1/4, 1/4),
		 w = c(.5, .5))

MW.nm8 <- norMix(name = "#8 Asym Bim",
		 mu   = c(0, 1.5),
		 sig2 = c(1, 1/9),
		 w = c(.75, .25))

MW.nm9 <- norMix(name = "#9 Trimodal",
		 mu   = c(-1.2, 1.2, 0),
		 sig2 = c(9/25, 9/25, 1/16),
		 w = c(.45, .45, .1))

MW.nm10 <- norMix(name = "#10 Claw",
		  mu = c(0, seq(-1,1, by = .5)),
		  sig2 = c(1, .01*rep(1,5)),
		  w = c(.5, .1*rep(1,5)))


MW.nm11 <- norMix(name = "#11 Doub Claw",
		  mu = c(-1, 1,seq(-1.5,1.5, by = .5)),
		  sig2 = c(4/9, 4/9, .0001*rep(1,7)),
		  w = c(.98*c(.5, .5), .02*rep(1,7)/7))


mu <- .5 + -2:2
sig2 <- c(.16, .04, .01, .0025, .000625)
w <- sqrt(sig2)
MW.nm12 <- norMix(name = "#12 Asym Claw",
		  mu   = c(0, mu),
		  sig2 = c(1, sig2),
		  w    = .5*c(1, w/sum(w)))


mu <- c(-1.5, -1, -.5)
sig2 <- rep(.0001, 3)
w <- rep(1/3,3)
MW.nm13 <- norMix(name = "#13 As Do Claw",
		  mu   = c(  -1,  1,   mu,  mu + 2),
		  sig2 = c(4/9, 4/9, sig2, 49*sig2),
		  w    = c(.46, .46, .01*w,  .07*w))


n <- 6
sig <- 2^(0:(1-n)) # =  2^{1-j} , j=1,2,..,n
il <- -3
##- ir <- 3
##- mu <- numeric(n)
##- for(i in 2:n) {
##-   mu[i] <- ir + 3 * sig[i] #-- mu_i  = ir_{i-1} + 3 s_i  = (ir_{i-1}+ir_i)/2
##-   ir    <- ir + 6 * sig[i] #-- ir_i  = -3 + sum_{j=1}^i (6 * s_j)
##- }
irv <- -3 + 6*cumsum(sig)
mu <- c(0, (irv[-1] + irv[-n])/2)
di <- irv[n] - il
MW.nm14 <- norMix(name = "#14 Smoo Comb",
		  mu   = 6 * ((mu - il) / di) - 3,
		  sig2 = (sig * 6 / di)^2,
		  w    = sig / sum(sig))


mu  <- c(0, 6, 12, 15.5, 16.5, 17.5)
sig <- c(1, 1,  1,  1/6,  1/6,  1/6)
il <- -3 ; ir <- 18
di <- ir - il
MW.nm15 <- norMix(name = "#15 Disc Comb",
		  mu   = 6 * ((mu - il) / di) - 3,
		  sig2 = (6 *sig/ di)^2,
		  w    = sig / sum(sig))

MW.nm16 <- norMix(name = "#16 Dist Bim",
		  mu   = c(-2.5, 2.5),
		  sig2 = c(1/36, 1/36),
		  w    = c(.5, .5))

## clean up:
rm(n, mu, sig, sig2, w, il,ir,irv,di)
