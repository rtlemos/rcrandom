# With this script you can confirm this versin of L'Ecuyer et al.'s MRG32k3a
# produces the same results as their C version.
# The output is supposed to emulate that of
# http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/test2RngStream.res
# 
# Author: Ricardo T Lemos
# Date: May 26, 2014
###############################################################################

cat("\nEmulating the tests performed by L'Ecuyer et al.
    in their C version of the RNG.\n")

sum3 <- sumi <- 0.0
if (exists("rcrandom.globalseed")) rm(rcrandom.globalseed)
if (exists("rcrandom.globalsubseed")) rm(rcrandom.globalsubseed)

# Three streams created with states that are 2^127 units apart
g1 <- rng() # If we do g1%advance.state(127,0), g1 reaches the state of g2
g2 <- rng()
g3 <- rng() # If we do g3%advance.state(-127,0), g3 reaches the state of g2

cat("Initial states of g1, g2 and g3:\n")
cat(g1$show(),"\n") # same as simply writing "g1",
cat(g2$show(),"\n")
cat(g3$show(),"\n")

summ <- g2$runif() + g3$runif()

cat("State of g1 after advancing by 2 ^ 5 + 3 = 35 steps:\n")
g1$advance.state(5,3)
cat(g1$show(),"\n")

cat("RandU01 (g1) =",g1$runif(),"\n")

#resetting the state of g1 to its original value (which is unique to stream g1)
g1$reset()
#cumbersome way to advance the state 35 times
for (i in 0:34) g1$advance.state(0, 1)
cat("State of g1 after reset and advancing 35 times by 1:\n")
cat(g1$show(),"\n")

cat("RandU01 (g1) =",g1$runif(),"\n")

g1$reset()
sumi <- sumi + sum(g1$rint(n = 35, lb = 1, ub = 10))
summ <- summ + sumi / 100
cat("State of g1 after reset and 35 calls to rint(1, 10):\n")
cat(g1$show(),"\n")
cat("   sum of 35 integers in [1,10]",sumi,"\n")

cat("RandU01 (g1) =",g1$runif(),"\n")

sumi <- 0.0
g1$reset()
g1$high.precision(TRUE)
sumi <- sumi + sum(g1$rint(n = 17, lb = 1, ub = 10))
cat("State of g1 after reset, IncreasedPrecis(T) and 17 calls to RandInt:\n")
cat(g1$show(),"\n")
g1$high.precision(FALSE)
g1$rint(1,10)
cat("State of g1 after IncreasedPrecis(F) and 1 calls to RandInt:\n")
cat(g1$show(),"\n")
sum3 <- sumi / 10.0

g1$reset()
g1$high.precision(TRUE)
sum3 <- sum3 + sum(g1$runif(n = 17))
cat("State of g1 after reset, IncreasedPrecis(T) and 17 calls to RandU01:\n")
cat(g1$show(),"\n")
g1$high.precision(FALSE)
z <- g1$runif()
cat("State of g1 after IncreasedPrecis(F) and 1 call to RandU01:\n")
cat(g1$show(),"\n")
summ <- summ + sum3 / 10.0

sum3 <- 0.0
sum3 <- sum3 + sum(g3$runif(n = 100))
cat("Sum of first 100 output values from stream g3:\n")
cat("   sum = ", sum3, "\n")
summ <- summ + sum3 / 10.0

cat("Resetting stream g3 to its initial seed.\n")
g3$reset()
cat("First 5 output values from stream g3:\n")
for (i in 1:5) cat(g3$runif(),"\n")
summ <- summ + g3$runif()

#Moving to the next substream; substreams are 2^74 states apart
cat("Resetting stream g3 to the next substream, 4 times.\n")
for (i in 1:4) g3$next.substream()
cat("First 5 output values from stream g3, fourth substream:\n")
for (i in 1:5) cat(g3$runif(),"\n")
summ <- summ + g3$runif()

cat("Resetting stream g2 to the beginning of substream.\n")
g2$reset()
g2$high.precision(TRUE)
sum3 <- 0.0
sum3 <- sum3 + sum(g2$runif(n = 100000))
summ <- summ + sum3 / 10000.0;
cat("Sum of 100000 values from stream g2 with double precision:",sum3,"\n")
g2$high.precision(FALSE)

g3$antithetic(TRUE)
sum3 <- 0.0
sum3 <- sum3 + sum(g3$runif(n = 100000))
summ <- summ + sum3 / 10000.0;
cat("Sum of 100000 antithetic values from stream g3:", sum3, "\n")

rm(g1, g2, g3)

cat("Setting global seed = { 1, 1, 1, 1, 1, 1 }\n")
rcrandom.globalseed    <- rep(1, 6)
rcrandom.globalsubseed <- rcrandom.globalseed

cat("Creating an array of 4 named streams and writing their full state:\n")
gar      <- vector("list",length = 4)
gar[[1]] <- rng(name = "Poisson")
gar[[2]] <- rng(name = "Laplace")
gar[[3]] <- rng(name = "Galois")
gar[[4]] <- rng(name = "Cantor")

for (i in 1:4) cat(gar[[i]]$show(complete = TRUE), "\n")

cat("Jumping stream Galois by 2^127 steps backward\n")
gar[[3]]$advance.state(-127, 0)
cat("The current state of the RNG Galois:\n")
cat(gar[[3]]$show(complete = TRUE), "\n")
gar[[3]]$next.substream()

for (i in 1:4) summ <- summ + gar[[i]]$runif()

cat("Final sum", summ, "\n")
