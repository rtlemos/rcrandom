% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r-mrg32k3a.R
\docType{class}
\name{rcmrg32k3a-class}
\alias{rcmrg32k3a}
\alias{rcmrg32k3a-class}
\title{rcmrg32k3a}
\description{
The RNG engine of this package.
}
\details{
This is an R+Fortran version of L'Ecuyer et al.'s MRG32k3a 
implementation in C, with a few extras for pseudo-random number generation
in random vectors. The class is just a shell for Fortran subroutines. It is
not exported because users of this package can employ the rcrng RC directly.

Citation

bibentry("article",
title = "An object-oriented random-number package with many long
streams and substreams",
author = "L'Ecuyer, P and Simard, R and Chen, EJ and Kelton, WD",
journal = "Operations Research",
year   = "2002",
address = "{ISI}:000179794700014",
volume   = "{ISBN} 3-900051-07-0",
number = "6",
pages    = "1073-1075",
mheader = "Paper describing the random number generator in C:")
}
\section{Fields}{

\describe{
\item{\code{algorithm.name}}{character. Name of algorithm}
}}
\section{Methods}{

\describe{
\item{\code{get.advseed(s)}}{Auxiliary function that provides the advanced seed}

\item{\code{get.advstate(ee, cc, cg)}}{Advance the state cg}

\item{\code{get.advsubseed(s, n.jumps = 1)}}{Auxiliary function that provides the advanced subseed}

\item{\code{get.nextsubstream(bg)}}{Obtain the state of the next substream}

\item{\code{get.seedcheck(seed)}}{Check if the seed is good}

\item{\code{get.u01(cg, anti, n = 1)}}{Obtain a uniform variate}

\item{\code{get.u01d(cg, anti, n = 1)}}{Obtain a high-precision uniform variate}

\item{\code{name(new.name = NULL)}}{Set name of algorithm object}
}}

