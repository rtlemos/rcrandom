#' rcmrg32k3a
#' 
#' The RNG engine of this package.
#' 
#' This is an R+Fortran version of L'Ecuyer et al.'s MRG32k3a 
#' implementation in C, with a few extras for pseudo-random number generation
#' in random vectors. The class is just a shell for Fortran subroutines. It is
#' not exported because users of this package can employ the rcrng RC directly.
#' 
#' Citation
#' 
#' bibentry("article",
#' title = "An object-oriented random-number package with many long
#' streams and substreams",
#' author = "L'Ecuyer, P and Simard, R and Chen, EJ and Kelton, WD",
#' journal = "Operations Research",
#' year   = "2002",
#' address = "{ISI}:000179794700014",
#' volume   = "{ISBN} 3-900051-07-0",
#' number = "6",
#' pages    = "1073-1075",
#' mheader = "Paper describing the random number generator in C:")
#'
#' @field algorithm.name character. Name of algorithm
#' 
#' @importFrom methods setRefClass
#' @export rcmrg32k3a
#' @exportClass rcmrg32k3a
#'
rcmrg32k3a <- setRefClass(
  Class = "rcmrg32k3a",
  fields = list(algorithm.name = "character"),
  methods = list(
    initialize = function(name = NULL){
      if (!missing(name)) .self$algorithm.name <- name
    },

    name = function(new.name = NULL){
      "Set name of algorithm object"

      if (!missing(new.name)) {
        .self$algorithm.name <- new.name
      } else {
        return(.self$algorithm.name)
      }
    },

    get.advseed = function(s){
      "Auxiliary function that provides the advanced seed"

      out <- .Fortran("get_advanced_seed", s = as.double(s) )
      return(out$s)
    },

    get.advsubseed = function(s, n.jumps = 1){
      "Auxiliary function that provides the advanced subseed"

      out <- .Fortran("get_advanced_subseed", s = as.double(s),
                      n_jumps = as.integer(n.jumps))
      return(out$s)
    },

    get.advstate = function(ee, cc, cg){
      "Advance the state cg"

      out <- .Fortran("get_advanced_state",
                      ee = as.integer(ee), cc = as.integer(cc),
                      Cg = as.double(cg))
      return(out$Cg)
    },

    get.seedcheck = function(seed){
      "Check if the seed is good"

      ok  <- 1
      out <- .Fortran("get_seed_check",
                      seed = as.double(seed),
                      ok = as.integer(ok))
      return(out$ok)
    },

    get.u01 = function(cg, anti, n = 1){
      "Obtain a uniform variate"

      s   <- ncol(cg)
      uu  <- rep(0.0, n * s)
      out <- .Fortran("U01", n = as.integer(n), s = as.integer(s),
                      Cg = as.double(cg), anti = as.logical(anti),
                      uu = as.double(uu))
      return(list(new.cg = matrix(nr = 6, nc = s, out$Cg), variate = out$uu))
    },

    get.u01d = function(cg, anti, n = 1){
      "Obtain a high-precision uniform variate"

      s   <- ncol(cg)
      uu  <- rep(0.0, n * s)
      out <- .Fortran("U01d", n = as.integer(n), s = as.integer(s),
                      Cg = as.double(cg), anti = as.logical(anti),
                      uu = as.double(uu))
      res <- list(new.cg = matrix(nr = 6, nc = s, out$Cg), variate = out$uu)
      return(res)
    },

    get.nextsubstream = function(bg){
      "Obtain the state of the next substream"

      out <- .Fortran("get_next_substream", Bg = as.double(bg))
      return(out$Bg)
    }
  )
)
