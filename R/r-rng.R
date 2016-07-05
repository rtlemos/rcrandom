#' rcrng
#' 
#' Pseudo-random number generator
#' 
#' This is the main reference class of this package. 
#' Use it to create pseudo-random number generators for random variables and
#' vectors. Currently it is assumed that the core algorithm
#' behind any rng object is Pierre L'Ecuyer's MRG32k3a. 
#' This assumption may be relaxed in future versions of this package.
#'
#' @field rng.name character. Name of the pseudo-random number generator object. 
#' @field n.substreams numeric. Number of substreams in the object. 
#' @field anti logical. Use anthitetical uniform variates (1-x) instead? 
#' @field inc.prec logical. Use variates with increased precision (53 bits)? 
#' @field resettable logical. Should the RNG be resettable?
#' @field cg matrix. State matrix 1.
#' @field bg matrix. State matrix 2. 
#' @field ig numeric. State matrix 3.
#' @field algorithm rcmrg32k3a. Algorithm that advances the state matrix. 
#' 
#' @importFrom methods new
#' @export rcrng
#' @exportClass rcrng
#'
#' @examples
#' h <- rcrng()
#' h$runif()
#' 
rcrng <- setRefClass(
  Class = "rcrng",
  fields = list(rng.name = "character", n.substreams = "numeric",
                anti = "logical", inc.prec = "logical", resettable = "logical",
                cg = "matrix", bg = "matrix", ig = "numeric",
                algorithm  = "rcmrg32k3a"),
  methods = list(
    
    initialize = function(n.substreams = 1, name = NULL, high.precision = FALSE,
                          antithetic = FALSE, resettable = TRUE, ...){
      "Creates a new random number generator.
      The actual algorithm may be created (default) or passed by reference."
      .self$anti <- antithetic
      .self$inc.prec <- high.precision
      .self$n.substreams <- n.substreams
      .self$resettable <- resettable
      
      if (!missing(name)) .self$rng.name <- name
      
      if (exists("rcrandom.globalalgorithm", where = .GlobalEnv)) {
        .self$algorithm <- rcrandom.globalalgorithm
      } else {
        .self$algorithm <- rcmrg32k3a()
      }
      
      #Generating seeds and subseeds if they have not been created yet
      if (!exists("rcrandom.globalseed", where = .GlobalEnv)) {
        gs <- rep(12345, 6)
        assign("rcrandom.globalseed", gs, envir = .GlobalEnv)
      }
      if (!exists("rcrandom.globalsubseed", where = .GlobalEnv)) {
        gs <- rep(12345, 6)
        assign("rcrandom.globalsubseed", gs, envir = .GlobalEnv)
      }
      
      # The following was introduced to check if substreams were created prior
      # to this call. If true, then the global seed must be advanced
      if (n.substreams == 1 &
          !identical(rcrandom.globalseed, rcrandom.globalsubseed)) {
        assign("rcrandom.globalseed",
               .self$algorithm$get.advseed(s = rcrandom.globalseed),
               envir = .GlobalEnv)
      }
      
      #initial seed for this RNG
      if (resettable) .self$ig <- rcrandom.globalseed
      
      #stream / substreams for this RNG
      if (n.substreams == 1) {
        .self$cg <- matrix(nrow = 6, ncol = 1, rcrandom.globalseed)
        #Advancing the seed & subseed for the next random variable
        advanced.seed <- .self$algorithm$get.advseed(s = rcrandom.globalseed)
        assign("rcrandom.globalseed", advanced.seed, envir = .GlobalEnv)
        assign("rcrandom.globalsubseed", advanced.seed, envir = .GlobalEnv)
      } else {
        .self$cg <- mapply(1:n.substreams, FUN = function(i){
          out <- rcrandom.globalsubseed
          assign("rcrandom.globalsubseed",
                 .self$algorithm$get.advsubseed(s = rcrandom.globalsubseed,
                                                n.jumps = 1),
                 envir = .GlobalEnv)
          return(out)
        })
      }
      #initial subseeds for this RNG
      if (resettable) .self$bg <- .self$cg
      
      callSuper(...)
    },
    
    runif = function(n = 1, n.ascolumn = TRUE){
      "Returns n pseudo-random numbers from the (0, 1) uniform distribution.
      In normal mode, the returned numbers have 32 bits of precision
      (i.e., they are multiples of 1 / (2 ^ 32 - 208)) and the
      state is advanced by one step per variate.
      In Increased Precision mode, returned numbers have 53 bits of precision,
      and the state is advanced by two steps per variate."
      
      if (.self$inc.prec) {
        res <- .self$algorithm$get.u01d(n = n, cg = .self$cg, anti = .self$anti)
      } else {
        res <- .self$algorithm$get.u01(n = n, cg = .self$cg, anti = .self$anti)
      }
      .self$cg <- res$new.cg
      if (.self$n.substreams == 1) {
        return(res$variate)
      } else if (n.ascolumn) {
        return(matrix(nrow = .self$n.substreams, ncol = n, res$variate,
                      byrow = TRUE))
      } else {
        return(matrix(nrow = n, ncol = .self$n.substreams, res$variate))
      }
    },
    
    rint = function(lb, ub, n = 1){
      "Returns a pseudo-random number from the discrete uniform distribution
      over the integers (i, i + 1, ... j). Makes one call to runif."
      
      u <- .self$runif(n = n)
      x <- as.integer(lb + as.integer(u * (ub - lb + 1)))
      return(x)
    },
    
    reset = function(){
      "Reinitializes the streams and substreams to their initial states, i.e.,
      if n = 1 (one stream only), cg and bg are set to ig;
      if n > 1 (multiple substreams), cg is set to bg."
      
      if (!.self$resettable) stop("Random number algorithm is not resettable.")
      if (.self$n.substreams == 1) {
        .self$cg[, 1] <- .self$bg[, 1] <- .self$ig
      } else {
        .self$cg <- .self$bg
      }
    },
    
    show = function(complete = FALSE){
      if (length(.self$rng.name) != 0) cat("The RNG", .self$rng.name, ":\n")
      if (complete & .self$resettable) {
        cat("ig = {", .self$ig, "}\n")
        for (i in 1:.self$n.substreams) {
          cat("bg[", i, "] = {", .self$bg[, i], "}\n")
        }
      }
      for (i in 1:.self$n.substreams) {
        cat("cg[", i,"] = {", .self$cg[, i], "}\n")
      }
    },
    
    seed = function(new.seed = NULL){
      "Provides current seed or sets a new seed, to produce the stream and
      substreams. Input value = vector with 6 integers."
      
      if (missing(new.seed)) {
        return(.self$ig) #returning the main seed
      }
      if (is(new.seed, "numeric") & length(new.seed) == 6) {
        myseed <- new.seed
      } else {
        stop("Unknown seed type.")
      }
      if (.self$algorithm$get.seedcheck(myseed) == 0) {
        for (i in 1:n.substreams) {
          .self$cg[, i] <- myseed
          myseed <- .self$algorithm$get.advsubseed(s = myseed, n.jumps = 1)
        }
        if (.self$resettable) {
          .self$ig <- myseed
          .self$bg <- .self$cg
        }
      } else {
        stop("Bad seed.")
      }
    },
    
    high.precision = function(value = NULL){
      "With value = TRUE, each call to the algorithm (direct or indirect) for
      this stream will return a uniform random number with 53 bits
      of resolution instead of 32 bits, and will advance the state
      of the stream by 2 steps instead of 1"
      if (missing(value)) {
        return(.self$inc.prec)
      } else {
        .self$inc.prec <- value
      }
    },
    
    antithetic = function(value = NULL){
      "If value = TRUE, the stream will start generating antithetic variates,
      i.e., 1-U instead of U.\n"
      
      if (missing(value)) {
        return(.self$anti)
      } else {
        .self$anti <- value
      }
    },
    
    next.substream = function(i = 1){
      "Moves a stream to the beginning of its next substream."
      
      if (!.self$resettable) {
        stop("This function is not available for this random number algorithm.")
      }
      .self$bg[, i] <- .self$cg[, i] <- .self$algorithm$get.nextsubstream(
        bg = .self$bg[, i])
    },
    
    advance.state = function(ee, cc, i = 1){
      "Advances the state by 2 ^ ee + cc units"
      
      .self$cg[, i] <- .self$algorithm$get.advstate(ee, cc, .self$cg[, i])
    }
    
  )
)
