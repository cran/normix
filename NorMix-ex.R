#### Examples of using MM's not-yet  "NorMix" Package :

## Load Martin Maechler's   'Normal Mixtures' package :
source("/u/maechler/R/MM/STATISTICS/NorMix.S")

## Testing:
nm <- NorMix(m= c(-1,1,2), sig2 = 3, w=1:3)
is.NorMix(nm)
nm
plot(dnm <- dNorMix(nm))

## Plot Method
plot(nm)
plot(NorMix(m= c(-1,1,2), sig2 = .8, w=1:3))
mult.fig(6,## main = "Normal Mix")
         main = expression("Normal Mixtures of 2 equals -- varying " * sigma))
for(s in 1.4^(-3:2))
    plot(NorMix(c(-1,1), sig2=s, long.name = TRUE),
         xlim=c(-6,6))


