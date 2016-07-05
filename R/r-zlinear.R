#' ZLinear
#'
#' Build piecewise linear log-probability density functions.
#'
#' This Reference Class is analogous to ZQuadratic.
#' Objects of this class provide approximations to the input
#' log-p.d.f. that are slightly more numerically stable,
#' but probably worse, than those of ZQuadratic.
#' To mitigate the worse fit, more evaluations should be
#' provided. For details about this classe's fields and
#' methods, see those in ZQuadratic.
#' Although this class is not virtual, it is not exported.
#'
#'
#' @field n numeric. 
#' @field x numeric. 
#' @field y numeric. 
#' @field beta0 numeric. 
#' @field beta1 numeric. 
#' @field dx numeric. 
#' @field normct numeric. 
#' @field pcdf numeric. 
#' @field fns list. 
#'
#' @importFrom methods new
#' @export ZLinear
#' @exportClass ZLinear
#'
#' @examples
#' x <- seq(-4, 4, 0.3)
#' y <- dnorm(x, log = TRUE)
#' zq <- ZLinear(x = x, y = y)
#' do.log <- TRUE #switch to FALSE to see p.d.f.
#' zq$plot.fit(do.log)
#' lines(seq(-4, 4, 0.1), dnorm(seq(-4, 4, 0.1),
#'  log = do.log), ty = "l", lwd = 1, lty = 2)
#'  
ZLinear <- setRefClass(
  Class = "ZLinear",
  fields = list(
    n = "numeric",
    x = "numeric",
    y = "numeric",
    beta0 = "numeric",
    beta1 = "numeric",
    dx = "numeric",
    normct = "numeric",
    pcdf = "numeric",
    fns = "list"),
  methods = list(
    initialize = function(x = NULL, y = NULL,
                          useFortran = FALSE){
      
      if (is.null(x) | is.null(y)) return()
      nx <- length(x)
      if (nx %% 2 != 1) {
        stop("Must use an odd number of evaluations.")
      }
      myf <- if (useFortran) get.ffns() else get.rfns()
      coeffs <- myf$betas.linear(nx, x, y)
      .self$x <- x #x-coordinates of evaluations
      .self$y <- y #these are the log-density evaluations
      .self$n <- nx
      .self$beta0 <- coeffs$b0
      .self$beta1 <- coeffs$b1
      .self$dx <- x[2:nx] - x[1:(nx - 1)]
      integrals <- myf$integrals.linear(coeffs$b0, coeffs$b1, 0, 1)
      .self$normct <- sum(integrals * .self$dx)
      .self$pcdf <- cumsum(integrals * .self$dx) / .self$normct
      .self$fns <- myf
    },
    
    pdf = function(z, do.log = FALSE){
      return(.self$fns$pdf.linear(
        z, .self$x, .self$dx, .self$beta0, .self$beta1,
        .self$normct, do.log, .self$fns$pdf.quadratic))
    },
    
    cdf = function(z){
      return(.self$fns$cdf.linear(
        z, .self$x, .self$dx, .self$beta0, .self$beta1,
        .self$pcdf, .self$normct))
    },
    
    invcdf = function(p){
      return(.self$fns$invcdf.linear(
        p, .self$x, .self$dx, .self$beta0, .self$beta1,
        .self$pcdf, .self$normct))
    },
    
    plot.fit = function(do.log = TRUE){
      if (do.log) {
        py <- .self$y
        mlab <- "piecewise linear approx. to log-pdf"
        yl <- "log-pdf"
      } else {
        py <- exp(.self$y)
        mlab <- "approximation to pdf"
        yl <- "pdf"
      }
      plot(.self$x, py, xlab = "z-coord", ylab = yl,
           main = mlab)
      for (i in 1:(.self$n - 1)) {
        fx <- seq(.self$x[i], .self$x[i + 1], length = 10)
        dd <- (fx - .self$x[i]) / .self$dx[i]
        fy <- .self$beta0[i] + .self$beta1[i] * dd
        ly <- if (do.log) fy else exp(fy)
        lines(fx, ly, col = i, lwd = 4)
      }
    },
    
    #-------------------------------------------------------
    # Private functions ------------------------------------
    #-------------------------------------------------------
    
    iget.coord = function(){
      return(.self$x)
    },
    
    iget.logd = function(){
      return(.self$y)
    },
    
    iget.multiply.constant = function(k){
      if (k > 0) {
        coords <- list(x = k * .self$x, y = .self$y)
      } else {
        coords <- list(x = k * rev(.self$x),
                       y = rev(.self$y))
      }
      return(coords)
    },
    
    is.operation.allowed = function(operation, argclass){
      switch( operation,
              "/" = FALSE,
              "*" = {any(argclass == c("numeric", "matrix",
                                       "Constant",
                                       "Uniform"))},
              "-" = ,
              "+" = {any(argclass == c("numeric",
                                       "matrix",
                                       "Constant",
                                       "Uniform",
                                       "Normal",
                                       "ZQuadratic",
                                       "ZLinear"))},
              FALSE #default
      )
    },
    iset.operate = function(operation, operand,
                            
                            operand.name, operand.side,
                            my.name){
      if (!.self$is.operation.allowed(operation,
                                      class(operand))) {
        stop("ZLinear: Invalid operation.")
      }
      my.specs <- list(n = length(.self$x), x = .self$x,
                       dx = .self$dx,
                       b0 = .self$beta0, b1 = .self$beta1,
                       b2 = rep(0,length(.self$beta0)),
                       normct = .self$normct)
      switch(
        class(operand),
        "ZLinear" = {
          op.specs <- list(
            n = length(operand$x),
            x = operand$x,
            dx = operand$dx,
            b0 = operand$beta0,
            b1 = operand$beta1,
            b2 = rep(0, length(operand$beta0)),
            normct = operand$normct)
          if (operand.side == "left") {
            coords <-
              .self$fns$convolution.zquadratic.zquadratic(
                l = op.specs,
                r = my.specs,
                operation = operation,
                find.cutpoints = .self$myf$find.cutpoints,
                integrals = .self$fns$integrals.linear,
                erf = .self$fns$erf,
                erfi = .self$fns$erfi)
          } else {
            coords <-
              .self$fns$convolution.zquadratic.zquadratic(
                l = op.specs,
                r = my.specs, operation = operation,
                find.cutpoints = .self$myf$find.cutpoints,
                integrals = .self$fns$integrals.linear,
                erf = .self$fns$erf,
                erfi = .self$fns$erfi)
          }
        },
        "ZQuadratic" = {
          op.specs <- list(n = length(operand$x),
                           x = operand$x,
                           dx = operand$dx,
                           b0 = operand$beta0,
                           b1 = operand$beta1,
                           b2 = operand$beta2,
                           normct = operand$normct)
          if (operand.side == "left") {
            coords <-
              .self$fns$convolution.zquadratic.zquadratic(
                l = op.specs, r = my.specs,
                operation = operation,
                find.cutpoints = .self$myf$find.cutpoints,
                integrals = .self$fns$integrals.linear,
                erf = .self$fns$erf,
                erfi = .self$fns$erfi)
          } else {
            coords <-
              .self$fns$convolution.zquadratic.zquadratic(
                l = op.specs,
                r = my.specs, operation = operation,
                find.cutpoints = .self$myf$find.cutpoints,
                integrals = .self$fns$integrals.linear,
                erf = .self$fns$erf,
                erfi = .self$fns$erfi)
          }
        },
        "Uniform" = {
          switch(
            operation,
            "-" = ,
            "+" = ,
            "*" = {
              unif.lb <- operand$param$lb$parameter(
                1, eval = TRUE)
              unif.ub <- operand$param$ub$parameter(
                1, eval = TRUE)
              unif.n <- operand$nr
              op.specs <- list(lb = unif.lb, ub = unif.ub,
                               n = unif.n )
              coords <-
                .self$fns$convolution.zquadratic.uniform(
                  zq = my.specs, unif = op.specs,
                  operation = operation,
                  find.cutpoints.zu =
                    .self$fns$find.cutpoints.zu,
                  integralsMult = .self$fns$integralsMult,
                  ExpIntegral = .self$fns$expint)
            },
            stop("ZLinear: Invalid operation, stage 2.")
          )
        },
        "Normal" = {
          mu <- operand$parameter(id = 1, eval = TRUE)
          tau <- operand$parameter(id = 2, eval = TRUE)
          nr <- operand$nr
          if (is(mu,"Constant")) {
            mu.val <- mu$parameter(id = 1, eval = TRUE)
          } else {
            stop("mean?")
          }
          if (is(tau,"Constant")) {
            tau.val <- tau$parameter(id = 1, eval = TRUE)
          } else {
            stop("var?")
          }
          op.specs <- list(m = mu.val, v = tau.val)
          if (operand.side == "left" & operation == "-") {
            # Normal - ZLinear
            tmp.coord <- .self$iget.multiply.constant(-1.0)
            tmp <- ZLinear(x = tmp.coord$x, y = tmp.coord$y)
            tmp.specs <- list(
              n = length(tmp$x),
              x = tmp$x,
              dx = tmp$dx,
              b0 = tmp$beta0,
              b1 = tmp$beta1,
              b2 = rep(0, length(tmp$beta0)),
              normct = tmp$normct)
            coords <-
              .self$fns$convolution.zquadratic.normal(
                zq = tmp.specs, normal = op.specs,
                operation = "+",
                integrals = .self$fns$integrals.linear,
                erf = .self$fns$erf,
                erfi = .self$fns$erfi)
          } else {
            coords <-
              .self$fns$convolution.zquadratic.normal(
                zq = my.specs, normal = op.specs,
                operation = operation,
                integrals = .self$fns$integrals.linear,
                erf = .self$fns$erf,
                erfi = .self$fns$erfi )
          }
        },
        stop("ZLinear: Invalid operation, stage 2") #default
      )
      return(coords)
    }
  )
)