# rcrandom: pseudo-random number generators for random variables and vectors
#
# This package lets you generate pseudo-random numbers using Pierre L'Ecuyer et
# al.'s MRG32K3a algorithm. The algorithm is coded in Fortran, for better
# performance. With this package you get independent streams for different
# random variables, as well as independent streams for the components of random
# vectors. To get started, type: 
# > vignette('rcrandom-vignette', package ='rcrandom')