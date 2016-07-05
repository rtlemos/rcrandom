#' opt
#'
#' Optimize parameters and approximate probability density functions of random
#' variables.
#'
#' Use objects of this class to optimize a p.d.f. and approximate it using
#' products of independent, transformed parameters. You need to initialize
#' an object of this class, call the constructor with a function and initial
#' guesstimates for the unknown parameters, and finally call the set.z method.
#' For convenience, you can do this in one step, when you initialize the object.
#' Then, use the pdf method to obtain fast evaluations of the approximate p.d.f.
#'
#' @field logpdf Log-pdf to be optimized
#' @field npar Number of unknown parameters
#' @field args Other arguments required by log-pdf
#' @field iteration Number of p.d.f. evaluations conducted
#' @field lower Lower bounds for unknown parameters
#' @field upper Upper bounds for unknown parameters
#' @field mode Mode of pdf (par. value that maximizes pdf)
#' @field fmode Maximum value of the log-pdf (value at mode)
#' @field lev0 Evaluation of the log-p.d.f. when pars = 0
#' @field covF Covariance factor at the mode
#' @field precR Precision factor at the mode
#' @field logdetPR Log determinant of par. cov. at the mode
#' @field optim.res Output from the optimization method
#' @field Z List of transformed model parameters
#'
#' @examples
#' myfun <- function(x) dgamma(x, shape = 1.3, scale = 1) *
#'                      dnorm(x, mean = 2, sd = 1)
#' myfun.integ <- integrate(
#'              f = myfun, lower = 0, upper = 10)$value
#' fn <- function(x) myfun(x) / myfun.integ
#' s <- opt(fn = fn, initial = 1, lower = 0)
#' x <- seq(0, 10, 0.1)
#' plot(x, fn(x), ty = "l", col = "grey", lwd = 4,
#'    xlab = "x", ylab = "pdf")
#' lines(x, s$pdf(x))
#' legend(0.6 * max(x), max(fn(x)),
#'        c("true pdf", "rcoptim pdf"),
#'        col = c("grey", "black"), lty = 1)
#'
#' @importFrom GenSA GenSA
#' @importFrom Matrix Diagonal
#' @importFrom Matrix chol
#' @importFrom Matrix diag
#' @importFrom Matrix solve
#' @importFrom Matrix t
#' @importFrom Matrix nearPD
#' @importFrom methods new
#' @export opt
#' @exportClass opt
#'
opt <- setRefClass(
  Class = "opt",
  fields = list(
    logpdf = "function",
    npar = "numeric",
    args = "list",
    iteration = "numeric",
    lower = "numeric", 
    upper = "numeric",
    z.lower = "numeric", 
    z.upper = "numeric",
    mode = "numeric", 
    fmode = "numeric",
    lev0 = "numeric", 
    covF = "Matrix",
    precR = "Matrix", 
    logdetPR = "numeric",
    optim.res = "list",
    Z = "list",
    verbose = "logical"
  ),
  methods = list(
    initialize = function(fn = NULL, initial = NULL,
                          is.logf = FALSE, lower = -Inf,
                          upper = Inf, max.it = 300,
                          args = NULL, zspecs = NULL,
                          skip.c = FALSE, skip.z = FALSE,
                          verbose = FALSE,
                          method = "GenSA",
                          optim.tol = 1e-4){
      "Lets you run optimization and pdf exploration in one go\n
      @param fn function p.d.f. (needn't integrate to 1) to be optimized\n
      @param initial numerical Initial values for the p.d.f. parameters\n
      @param is.logf logical Is fn already log-transformed?\n
      @param lower numerical Lower bounds for the parameters\n
      @param upper numerical Upper bounds for the parameters\n
      @param args list Additional input arguments for the function\n
      @param zspecs list Specifications used to transform the parameters\n
      @param skip.z logical Should parameter transformation be skipped? \n
      @param method character Optimization algorithm \n
      @param optim.tol numerical Optimization relative tolerance.\n
      "

      if (skip.c | is.null(fn) | is.null(initial)) {
        return()
      } else if (is.null(fn) | is.null(initial)) {
        stop("To construct this object, input function and
             initial parameter estimates must be provided.")
      }
      .self$construct(
        fn = fn, is.logf = is.logf, initial = initial,
        lower = lower, upper = upper, max.it = max.it,
        args = args, zspecs = zspecs, skip.z = skip.z,
        verbose = verbose, method = method,
        optim.tol = optim.tol)
      },

    construct = function(
      fn, initial, is.logf = FALSE, lower = -Inf,
      upper = Inf, max.it = 300, args = NULL, zspecs = NULL,
      skip.z = FALSE, verbose = FALSE,
      method = "Nelder-Mead", optim.tol = 1e-4){
      "Optimization method that lets you find the mode of a p.d.f. \n
      @param fn function p.d.f. (needn't integrate to 1) to be optimized\n
      @param is.logf logical Is fn already log-transformed?\n
      @param initial numerical Initial values for the p.d.f. parameters\n
      @param lower numerical Lower bounds for the parameters\n
      @param upper numerical Upper bounds for the parameters\n
      @param args list Additional input arguments for the function\n
      @param zspecs list Specifications used to transform the parameters\n
      @param skip.z logical Should parameter transformation be skipped?
      @param method character Optimization algorithm
      @param optim.tol numerical Optimization relative tolerance
      "

      .self$npar <- length(initial)
      .self$iteration <- 0
      .self$verbose <- verbose
      n <- (2 ^ .self$npar + 1) * max.it
      if (is.list(args)) {
        .self$args <- args
      } else {
        .self$args <- vector("list")
      }
      if (length(lower) == .self$npar) {
        .self$lower <- lower
      } else {
        .self$lower <- rep(lower, .self$npar)
      }
      if (length(upper) == .self$npar) {
        .self$upper <- upper
      } else {
        .self$upper <- rep(upper, .self$npar)
      }
      if (any(names(formals(fn)) == "args")) {
        ev0 <- fn(rep(0, .self$npar), args = .self$args)
      } else {
        ev0 <- fn(rep(0, .self$npar))
      }
      .self$lev0 <- if (is.logf) ev0 else log(ev0)
      .self$logpdf <- function(x, update.counter = TRUE){
        if (identical(x, rep(0, .self$npar))) {
          lev <- .self$lev0
        } else {
          if (any(names(formals(fn)) == "args")) {
            ev <- fn(x, args = .self$args)
          } else {
            ev <- fn(x)
          }
          lev <- if (is.logf) ev else log(ev)
          if (update.counter) {
            .self$iteration <- .self$iteration + 1
            if (.self$verbose) {
              if (.self$iteration %% 50 == 0) {
                cat(paste0("*", formatC(.self$iteration, width = 8)), "\n")
              } else {
                cat("*")
              }
            }
          }
        }
        return(lev)
      }
      if (all(is.null(initial))) {
        ini <- rep(0, .self$npar)
      } else {
        ini <- initial
      }

      #Finding the optimum
      get.scaled.evaluation <- function(x, update.counter){
        ev <- .self$logpdf(x = x, update.counter = update.counter)
        return(-1 * ev)
      }
      .self$optim.res <- switch(
        method,
        "GenSA" = GenSA::GenSA(
          par = ini, fn = get.scaled.evaluation,
          lower = .self$lower, upper = .self$upper,
          control = list(max.call = max.it),
          update.counter = TRUE),
        "Nelder-Mead" = ,
        "BFGS" = ,
        "CG" = optim(par = ini,
                     fn = get.scaled.evaluation,
                     method = method,
                     control = list(reltol = optim.tol),
                     update.counter = TRUE),
        "L-BFGS-B" = optim(
          par = ini,fn = get.scaled.evaluation,
          method = method, lower = .self$lower,
          upper = .self$upper,
          control = list(reltol = optim.tol),
          update.counter = TRUE)
      )
      .self$mode <- .self$optim.res$par
      .self$fmode <- -1 * .self$optim.res$value
      if (!skip.z) {
        .self$set.z(zspecs = zspecs)
        if (.self$verbose) {
          cat("\n")
          .self$show()
        }
      }
    },

    set.z = function(zspecs = NULL){
      "Provides transformed input parameters based on the optimization of
      the input function.\n
      @param zspecs list list of specifications, e.g.\n
      zspecs <- list(nside = 10, logtol = -5, eps = 1e-4,\n
      useFortran = TRUE, useQuadratic = TRUE, useEigen = TRUE)
      "

      ub <- .self$upper - 1e-8
      lb <- .self$lower + 1e-8
      if (is.null(zspecs$nside)) {
        nside <- 10
      } else {
        nside <- zspecs$nside
      }
      if (is.null(zspecs$logtol)) {
        logtol <- -5
      } else {
        logtol <- zspecs$logtol
      }
      eps <- if (is.null(zspecs$eps)) 1e-4 else zspecs$eps
      if (is.null(zspecs$useEigen)) {
        useEigen <- TRUE
      } else {
        useEigen <- zspecs$useEigen
      }
      if (is.null(zspecs$useFortran)) {
        useFortran <- TRUE
      } else {
        useFortran <- zspecs$useFortran
      }
      if (is.null(zspecs$useQuadratic)) {
        useQuadratic <- TRUE
      } else {
        useQuadratic <- zspecs$useQuadratic
      }

      x1 <- mapply(1:.self$npar, FUN = function(i){
        x <- .self$mode
        x[i] <- max(lb[i], .self$mode[i] - eps)
        return(x)
      })
      x2 <- mapply(1:.self$npar, FUN = function(i){
        x <- .self$mode
        x[i] <- min(ub[i], .self$mode[i] + eps)
        return(x)
      })
      x3 <- if (.self$npar > 1) {
        matrix(
          nrow = .self$npar,
          unlist(mapply(1:(.self$npar - 1),
                        FUN = function(i){
                          mapply((i + 1):.self$npar,
                                 FUN = function(j){
                                   x <- .self$mode
                                   x[i] <- min(ub[i], .self$mode[i] + eps)
                                   x[j] <- min(ub[j], .self$mode[j] + eps)
                                   return(x)
                                 })
                        })))
      }
      x <- t(cbind(x1, x2, x3))

      yM <- .self$fmode
      yL <- mapply(1:.self$npar, FUN = function(i){
        x <- .self$mode
        x[i] <- max(lb[i], .self$mode[i] - eps)
        ev <- .self$logpdf(x = x, update.counter = TRUE)
        return(ev)
      })
      yU <- mapply(1:.self$npar, FUN = function(i){
        x <- .self$mode
        x[i] <- min(ub[i], .self$mode[i] + eps)
        ev <- .self$logpdf(x = x, update.counter = TRUE)
        return(ev)
      })
      if (.self$npar > 1) {
        yC <- lapply(1:(.self$npar - 1), FUN = function(i){
          k <- i + 1
          vec <- mapply(k:.self$npar, FUN = function(j){
            x <- .self$mode
            x[i] <- min(ub[i], .self$mode[i] + eps)
            x[j] <- min(ub[j], .self$mode[j] + eps)
            ev <- .self$logpdf(x = x, update.counter = TRUE)
            return(ev)
          })
          return(vec)
        })
      } else {
        yC <- NULL
      }

      lev <- c(yL, yU, unlist(yC))
      if (any(is.na(lev))) {
        stop("Some evaluations are NaNs. Check parameter bounds.")
      }
      if (any(is.infinite(lev))) {
        stop("Some evaluations are Inf. Check parameter bounds.")
      }

      prec <- mapply(1:.self$npar, FUN = function(i){
        vec <- mapply(1:.self$npar, FUN = function(j){
          if (i == j) {
            b <- (0.5 * (yL[i] + yU[i]) - yM) / eps ^ 2
            out <- -2 * b
          } else {
            yij <- yC[[min(i, j)]][max(i, j) - min(i,j)]
            b <- (yij + yM - (yU[i] + yU[j])) / eps ^ 2
            out <- -b
          }
          return(out)
        })
        return(vec)
      })

      if (.self$npar > 1) {
        p.aux <- Matrix::nearPD(prec)
      } else {
        p.aux <- list(mat = prec)
      }
      if (useEigen) {
        e <- eigen(p.aux$mat)
        .self$precR <- Matrix::Diagonal(
          n = .self$npar,
          x = sqrt(e$values)) %*% t(e$vectors)
        .self$logdetPR <- 0.5 * sum(log(e$values))
        .self$covF <- e$vectors %*% Matrix::Diagonal(
          n = .self$npar, x = 1 / sqrt(e$values))
      } else {
        .self$precR <- Matrix::chol(p.aux$mat)
        .self$logdetPR <- sum(log(
          Matrix::diag(.self$precR)))
        .self$covF <- Matrix::solve(.self$precR)
      }

      #Finding the bounds of the transformed variables
      xbnd.transform <- function(xbnd){
        out <- mapply(1:.self$npar, FUN = function(j){
          vec <- rep(0, .self$npar)
          vec[j] <- xbnd[j]
          res <- .self$get.x2z(vec)
          return(res)
        })
        return(out)
      }
      zbnd <- lapply(1:2, FUN = function(i){
        xbnd <- if (i == 1) .self$lower else .self$upper
        res <- xbnd.transform(xbnd)
      })
      if (.self$npar == 1) {
        .self$z.lower <- min(zbnd[[1]], zbnd[[2]])
        .self$z.upper <- max(zbnd[[1]], zbnd[[2]])
      } else {
        .self$z.lower <-
          mapply(1:.self$npar, FUN = function(i){
            vec <- mapply(1:.self$npar, FUN = function(j){
              min(zbnd[[1]][i,j], zbnd[[2]][i,j])
            })
            return(max(vec))
          })
        .self$z.upper <-
          mapply(1:.self$npar, FUN = function(i){
            vec <- mapply(1:.self$npar, FUN = function(j){
              max(zbnd[[1]][i,j], zbnd[[2]][i,j])
            })
            return(min(vec))
          })
      }

      #creating transformed variables
      zsamples <- .self$get.x2z(t(x))
      search.fun <- function(z, i, j){
        vec <- rep(0, .self$npar)
        vec[i] <- if (j == 1) -exp(z) else exp(z)
        x <- .self$get.z2x(vec)
        lev <- .self$logpdf(x, update.counter = TRUE)
        out <- abs(lev - .self$fmode - logtol)
        return(out)
      }
      .self$Z <- lapply(1:.self$npar, FUN = function(i){

        new.z <- mapply(1:2, FUN = function(j){
          opt <- GenSA::GenSA(
            par = -5, fn = search.fun,
            i = i, j = j, lower = -100, upper = 100,
            control = list(threshold.stop = 0.001, max.time = 2, 
                           trace.mat = FALSE))
          res <- if (j == 1) -exp(opt$par) else exp(opt$par)
          return(res)
        })

        # location of points for exact density evaluation
        if (nside < 10) {
          decay <- -log(0.01) / (nside - 1)
          ptsL <- exp(-decay * (0:(nside - 1)))
          ptsR <- exp(-decay * ((nside - 1):0))
        } else {
          ptsR <- seq(1 / nside, 1, length = nside)
          ptsL <- rev(ptsR)
        }
        final.z <- c(new.z[1] * ptsL, 0, new.z[2] * ptsR)
        vec <- rep(0, .self$npar)
        final.ld <- mapply(final.z, FUN = function(myz){
          vec[i] <- myz
          x <- .self$get.z2x(z = vec)
          ld <- .self$logpdf(x = x, update.counter = TRUE)
          return(ld)
        })

        #using ZQuadratic or Zlinear
        if (useQuadratic) {
          out <- ZQuadratic(x = final.z,
                            y = final.ld,
                            useFortran = useFortran)
        } else {
          out <- ZLinear(x = final.z,
                         y = final.ld,
                         useFortran = useFortran)
        }
        return(out)
      })
      if (.self$verbose) cat("\n")
      },

    get.x2z = function(x){
      "Provides transformed parameters based on
      original input parameters.\n
      @param x numeric, matrix, Matrix X-coordinates to be
      transformed\n
      @result z numeric, matrix, Matrix Transformed
      z-coordinates
      "

      z <- switch(
        class(x),
        "matrix" = t(as.matrix(.self$precR) %*%
                       (x - .self$mode)),
        "Matrix" = Matrix::t(.self$precR %*%
                               (x - .self$mode)),
        "numeric" = as.numeric(.self$precR %*%
                                 (x - .self$mode)),
        NA
      )
      return(z)
    },

    get.z2x = function(z){
      "Provides parameters with their original units based
      on transformed parameters.\n
      @param z numeric, matrix, Matrix Z-coordinates to be
      transformed\n
      @result x numeric, matrix, Matrix Transformed
      x-coordinates
      "

      x <- switch(
        class(z),
        "matrix" = as.matrix(.self$covF) %*% z + .self$mode,
        "Matrix" = .self$covF %*% z + .self$mode,
        "numeric" = as.numeric(.self$covF %*% z) +
          .self$mode,
        NA
      )
      return(x)
    },

    get.zaxis = function(id = 1){
      "Coordinates of transformed parameters,
      at all evaluation points.\n
      @param id numeric Identification number of Z-axis\n
      @result axis matrix X-coordinates of required Z-axis
      "

      if (length(.self$Z) == 0) .self$set.z()
      zcoord <- .self$Z[[id]]$x
      mat <- matrix(nrow = .self$npar,
                    ncol = length(zcoord), 0)
      mat[id, ] <- zcoord
      axis <- t(.self$get.z2x(mat))
      return(axis)
    },

    iget.reshape = function(f, x = NULL){
      "Internal function that reshapes the input provided
      to pdf, cdf and invcdf"

      if (is.null(x)) {
        out <- f
      } else {
        if (is.numeric(x)) {
          out <- f(matrix(nrow = .self$npar, x))
        } else if (ncol(x) == .self$npar &
                   nrow(x) != .self$npar) {
          out <- f(t(x))
        }
      }
      return(out)
    },

    get.extremes = function(id){
      "Provides the bounds of the region where a given input
      parameter has non-zero density"

      x1 <- mapply(1:.self$npar, FUN = function(i){
        vec <- rep(0, .self$npar)
        vec[i] <- .self$Z[[i]]$x[1]
        x <- .self$get.z2x(vec)
        return(x[id])
      })
      x2 <- mapply(1:.self$npar, FUN = function(i){
        vec <- rep(0, .self$npar)
        vec[i] <- .self$Z[[i]]$x[.self$Z[[i]]$n]
        x <- .self$get.z2x(vec)
        return(x[id])
      })
      lb <- min(x1, x2)
      ub <- max(x1, x2)
      return(c(lb, ub))
    },

    pdf = function(x = NULL, do.log = FALSE){
      "Approximate p.d.f. evaluation.\n
      @param x numeric,matrix,NULL X-coords of points where pdf is evaluated\n
      @param do.log logical Provide log-transformed p.d.f. instead?\n
      @result myf numeric,matrix,function Either evaluations or R function
      "

      if (length(.self$Z) == 0) .self$set.z()
      myf <- function(x){
        z <- .self$get.x2z(x)
        raw <- matrix(
          ncol = .self$npar,
          mapply(1:.self$npar, FUN = function(j) {
            .self$Z[[j]]$pdf(z = z[,j], do.log = do.log)
          }))
        if (do.log) {
          eval <- apply(raw, MARGIN = 1, FUN = sum)
          out <- eval + .self$logdetPR
        } else {
          eval <- apply(raw, MARGIN = 1, FUN = prod)
          out <- eval * exp(.self$logdetPR)
        }
        return(out)
      }
      out <- .self$iget.reshape(f = myf, x = x)
      return(out)
    },

    cdf = function(x = NULL, ...){
      "Provides either cumulative density function or its evaluation given x"

      myf <- function(x){
        z <- .self$get.x2z(x)
        raw <- matrix(
          ncol = .self$npar,
          mapply(1:.self$npar, FUN = function(i){
            .self$Z[[i]]$cdf(z[, i])
          }))
        if (is(raw,'numeric')) {
          eval <- raw
        } else {
          eval <- apply(raw, MARGIN = 1, FUN = prod)
        }
        return(eval)
      }
      out <- .self$iget.reshape(f = myf, x = x)
      return(out)
    },

    invcdf = function(p = NULL, ...){
      "Provides either inverse of the cumulative density
      function or its evaluation, for a given p."

      if (.self$npar > 1) {
        stop("Function unavailable for multivariate objs.")
      }
      if (!is.null(p)) if (any(p < 0) | any(p > 1)) {
        stop("Probability must be between 0 and 1")
      }
      myf <- function(p){
        z <- .self$Z[[1]]$invcdf(p)
        out <- .self$get.z2x(z)
        return(out)
      }
      out <- .self$iget.reshape(f = myf, x = p)
      return(out)
    },

    #Method that replaces default
    show = function(){
      "Summary of object"

      if (.self$npar == 1) {
        cat(paste0('Univariate Continuous (mode_x = ', .self$mode,
            ', mode_var = ', .self$covF[1,1] ^ 2,')\n'))
      } else {
        cat(paste0("Multivariate Continuous ",
                   "with coordinates at mode\n"))
        print(.self$mode)
        cat("and covariance at mode\n")
        print(crossprod(.self$covF))
      }
    }
  )
)