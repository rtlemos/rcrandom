## ---- warning = FALSE, message = FALSE-----------------------------------
library(rcrandom)
# ------------------------------------------------------------------------------
# Simplest case: one random number generator (RNG) with one stream
# ------------------------------------------------------------------------------
g <- rcrng()
# Generating 4 uniform variates with this RNG
g$runif(n = 4)
# Example of reproducibility: resetting stream and regenerating 4 variates
g$reset()
g$runif(n = 4)
#Generating integer uniform variates
g$rint(n = 12, lb = -3, ub = 3)
# A new generator has a different initial seed (2^127 units away from 1st one); 
# thus, random numbers differ from the previous generator's.
k <- rcrng()
k$runif(n = 4)
# A 3rd generator, for a random vector of size 2, generates 2 variates per call
h <- rcrng(2)
h$runif()
# Resetting this RNG and requesting 3 pairs of variates
h$reset()
h$runif(n = 3)
# ------------------------------------------------------------------------------
# Special feature of R Reference Classes: using 2 pointers to the same generator
# ------------------------------------------------------------------------------
pointer1 <- rcrng()
pointer2 <- pointer1
pointer1 #this shows the state of the RNG
pointer2
pointer1$runif(3) #generating variates with pointer1
pointer1
pointer2 #having used pointer1 to get variates also changed state of pointer2
pointer2$runif(1) #using pointer2 also changes pointer1
pointer1
pointer2

## ---- warning = FALSE, message = FALSE-----------------------------------
# ------------------------------------------------------------------------------
# Memory saving strategies
# ------------------------------------------------------------------------------
# Sharing the same algorithm across RNGs and disabling the "reset button" 
rcrng.globalalgorithm  <- rcmrg32k3a(name = "my.algo")
g1 <- rcrng(resettable = FALSE) #this RNG can't be reset
g2 <- rcrng()
g1$algorithm$name('new.name') #not a common procedure
rcrng.globalalgorithm$name() #name was changed
g1$runif()
g2$runif()
# g1$reset() #would throw an error
g2$reset()
g1$runif() #g1 was not reset
g2$runif() #g2 was reset
# ------------------------------------------------------------------------------
# Linked pseudo-random number generators
# ------------------------------------------------------------------------------
# a) Creating 3 RNGs with the same seed but different modes
h1 <- rcrng()
h2 <- rcrng()
h2$seed(h1$seed()) # one way of setting the seed of one RNG equal to another
h3 <- h1$copy()    # this is a more convenient way to do the same thing
h2$antithetic(TRUE)
h3$high.precision(TRUE)
h1
h2
h3
h1$runif() # u
h2$runif() # 1 - u 
h3$runif() # a different variate
h1
h2 #this stream is still aligned with h1's
h3 #this stream moves twice as fast as the previous ones
# b) Creating lagged RNGs
k1 <- rcrng(1)
k2 <- k1$copy()
k3 <- k1$copy()
k2$advance.state(ee = 0, cc = 1) #advancing the state by 2^ee+cc = 1 
k3$advance.state(ee = 1, cc = 0) #advancing by 2
k1
k2
k3
k1$runif(3)
k2$runif(3)
k3$runif(3)
# Confirming that the substream responsible for the 2nd element in the random 
# vector is 2^76 (approximately 7.6*10^22) units apart from the 1st element
j <- rcrng(2)
j$runif(n = 3)
j$next.substream(1)
j$reset()
j$runif(n = 3)

## ---- warning = FALSE, message = FALSE-----------------------------------
library("snowfall")
sfInit(parallel = TRUE, cpus = 2) 
sfLibrary(rcrandom) #slave processes load package 
# ------------------------------------------------------------------------------
# Reproducible variates generated in parallel
# ------------------------------------------------------------------------------
# Creating a local RNG
g <- rcrng(3)
sfExport("g")
# Generating 5 triplets of variates locally
g$runif(5)
#generating the same 5 triplets of variates remotely (on the 2 slave processes)
sfLapply(1:2, function(i) g$runif(5))
# ------------------------------------------------------------------------------
# Lagged variates (lags = 0, 1, 2) generated in parallel
# ------------------------------------------------------------------------------
k <- rcrng(3)
sfExport("k")
ok <- sfLapply(0:1, function(l) lapply(1:3, function(i){
  k$advance.state(ee = 0, cc = l * (i - 1), i = i)
}))
sfLapply(1:2, function(i) k$runif(5))
# ------------------------------------------------------------------------------
# Uncorrelated variates from 2 slave processes
# ------------------------------------------------------------------------------
h <- rcrng()
sfExport("h")
ok <- sfLapply(1:2, function(i) {if (i == 2) h$next.substream(); return()} )
sfLapply(1:2, function(i) h$runif(5))
# ------------------------------------------------------------------------------
# Terminating slave processes
# ------------------------------------------------------------------------------
sfStop()

