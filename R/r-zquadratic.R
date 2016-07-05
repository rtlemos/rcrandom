#' ZQuadratic
#'
#' Build piecewise quadratic log-probability density functions
#'
#' This is an auxiliary Reference Class in this package.
#' Objects of this reference class are created by opt objects, to represent 
#' the p.d.f. of transformed, orthogonal parameters. One object is created per 
#' p.d.f. parameter. Users of the package needn't interact with these
#' objects directly; instead, they should call opt object methods. 
#' PQFA = piecewise quadratic approximation to log-probability density function.
#'
#' @field n numeric No. points of exact p.d.f. evaluations
#' @field x numeric Univariate z-coordinates of the n points
#' @field y numeric Exact p.d.f. evaluations of the n points
#' @field beta0 numeric Intercepts of PQFA
#' @field beta1 numeric Slopes of PQFA
#' @field beta2 numeric Coeffs. of quadratic terms in PQFA
#' @field dx numeric Distances between adjac. z-coordinates
#' @field normct numeric Normalizing constant of PQFA
#' @field pcdf numeric Partial cdf associated with PQFA
#' @field fns list Auxiliary functions (in R or Fortran)
#'
#' @importFrom methods new
#' @export ZQuadratic
#' @exportClass ZQuadratic
#'
#' @examples
#' x <- c(-4:4)
#' y <- dnorm(x, log = TRUE)
#' zq <- ZQuadratic(x = x, y = y)
#' do.log <- TRUE #switch to FALSE to see p.d.f.
#' zq$plot.fit(do.log)
#' lines(seq(-4, 4, 0.1), dnorm(seq(-4, 4, 0.1),
#'       log = do.log), ty = "l", lwd = 1, lty = 2)
#'
ZQuadratic <- setRefClass(
  Class = "ZQuadratic",
  fields = list(
    n = "numeric",
    x = "numeric",
    y = "numeric",
    beta0 = "numeric",
    beta1 = "numeric",
    beta2 = "numeric",
    dx = "numeric",
    normct = "numeric",
    pcdf = "numeric",
    fns = "list"),
  methods = list(
    initialize = function(x = NULL, y = NULL, useFortran = TRUE){
      "Create a new ZQuadratic object.\n
      @param x numeric Vector of n z-coordinates, n odd, x[(n-1)/2 + 1] = 0 \n
      @param y numeric Vector of n exact log-pdf evaluations, made at x \n
      @param useFortran Should Fortran fns. be used instead of R fns.?
      "
      
      if (is.null(x) | is.null(y)) return()
      nx <- length(x)
      if (nx %% 2 != 1) {
        stop("Must provide an odd number of evaluations.")
      }
      myf <- if (useFortran) get.ffns() else get.rfns()
      coeffs <- myf$betas.quadratic(nx, x, y)
      .self$x <- x #x-coordinates of evaluations
      .self$y <- y #these are the log-density evaluations
      .self$n <- nx
      .self$beta0 <- coeffs$b0
      .self$beta1 <- coeffs$b1
      .self$beta2 <- coeffs$b2
      .self$dx <- x[2:nx] - x[1:(nx - 1)]
      integrals <- myf$integrals.quadratic(
        coeffs$b0, coeffs$b1, coeffs$b2, 0, 1,
        myf$erf, myf$erfi)
      .self$normct <- sum(integrals * .self$dx)
      .self$pcdf <- cumsum(integrals * .self$dx) / .self$normct
      .self$fns <- myf
      #.self$operations.classes <- list("/"=NULL,
      # "*"=c("numeric","matrix","Constant","Uniform"),
      # "-"=c("numeric","matrix","Constant","Uniform",
      #       "Normal", "ZQuadratic"),
      # "+"=c("numeric","matrix","Constant","Uniform",
      #       "Normal", "ZQuadratic"))
    },
    
    pdf = function(z, do.log = FALSE){
      "Probability density function of the transformed z-parameter\n
      @param z numeric Input (transformed) parameter for the approx. pdf\n
      @param do.log logical Log-transform the approx. pdf evaluation?\n
      @result Evaluation(s) of the approximate pdf
      "
      
      return(.self$fns$pdf.quadratic(
        z, .self$x, .self$dx, .self$beta0,
        .self$beta1, .self$beta2,
        .self$normct, do.log))
    },
    
    cdf = function(z){
      "Cumulative distribution function of the transformed z-parameter.\n
      @param z numeric Input (transformed) parameter for the approx. pdf\n
      @result Evaluation(s) of the approximate c.d.f.
      "
      
      return(.self$fns$cdf.quadratic(
        z, .self$x, .self$dx, .self$beta0,
        .self$beta1, .self$beta2, .self$pcdf,
        .self$normct, .self$fns$erf, .self$fns$erfi))
    },
    
    invcdf = function(p){
      "Inverse c.d.f. of the input probability.\n
      @param p numeric A probability\n
      @result Evaluation(s) of the approximate inverse c.d.f.
      "
      
      return(.self$fns$invcdf.quadratic(
        p, .self$x, .self$dx, .self$beta0, .self$beta1,
        .self$beta2, .self$pcdf, .self$normct,
        .self$fns$erf, .self$fns$erfi,
        .self$fns$inv.erf, .self$fns$inv.erfi))
    },
    
    plot.fit = function(do.log = TRUE){
      "Plot the piecewise approximation to the p.d.f.
      @param do.log logical Display log-p.d.f.?
      "
      
      if (do.log) {
        py <- .self$y
        yl <- "log-pdf"
      } else {
        py <- exp(.self$y)
        yl <- "pdf"
      }
      plot(.self$x, py, xlab = "z-coord", ylab = yl)
      for (i in 1:(.self$n - 1)) {
        fx <- seq(.self$x[i], .self$x[i + 1], length = 10)
        dd <- (fx - .self$x[i]) / .self$dx[i]
        fy <- .self$beta0[i] + .self$beta1[i] * dd +
          .self$beta2[i] * dd ^ 2
        ly <- if (do.log) fy else exp(fy)
        lines(fx, ly, col = i, lwd = 4)
      }
    },
    
    #-------------------------------------------------------
    # Private functions ------------------------------------
    #-------------------------------------------------------
    
    iget.coord = function(){
      "Returns transformed z-coordinates"
      
      return(.self$x)
    },
    
    iget.logd = function(){
      "Returns exact p.d.f. evaluations"
      
      return(.self$y)
    },
    
    iget.multiply.constant = function(k){
      "Multiplies random variable (Z-parameter) by constant"
      
      if (k > 0) {
        coords <- list(x = k * .self$x, y = .self$y)
      } else {
        coords <- list(x = k * rev(.self$x),
                       y = rev(.self$y))
      }
      return(coords)
    },
    
    is.operation.allowed = function(operation, argclass){
      "Checks if arithmetic operation is allowed"
      
      switch(
        operation,
        "/" = FALSE,
        "*" = {any(argclass == c("numeric", "matrix", "Constant", "Uniform"))},
        "-" = ,
        "+" = {any(argclass == c("numeric", "matrix", "Constant", "Uniform",
                                 "Normal", "ZQuadratic"))},
        FALSE #default
      )
    },
    
    iset.operate = function(operation, operand,
                            operand.name, operand.side,
                            my.name){
      "Performs arithmetic operation on random variable (Z-parameter)"
      
      if (!.self$is.operation.allowed(operation, class(operand))) {
        stop("ZQuadratic: Invalid operation.")
      }
      my.specs <- list(n = length(.self$x), x = .self$x,
                       dx = .self$dx, b0 = .self$beta0,
                       b1 = .self$beta1, b2 = .self$beta2,
                       normct = .self$normct)
      switch(
        class(operand),
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
                integrals = .self$fns$integrals.quadratic,
                erf = .self$fns$erf, erfi = .self$fns$erfi)
          } else {
            coords <-
              .self$fns$convolution.zquadratic.zquadratic(
                l = op.specs, r = my.specs,
                operation = operation,
                find.cutpoints = .self$myf$find.cutpoints,
                integrals = .self$fns$integrals.quadratic,
                erf = .self$fns$erf, erfi = .self$fns$erfi)
          }
        },
        "Uniform" = {
          switch(
            operation,
            "-" = ,
            "+" = ,
            "*" = {
              unif.lb <- operand$param$lb$parameter(1, eval = TRUE)
              unif.ub <- operand$param$ub$parameter(1, eval = TRUE)
              unif.n <- operand$nr
              op.specs <- list(lb = unif.lb, ub = unif.ub, n = unif.n )
              coords <-
                .self$fns$convolution.zquadratic.uniform(
                  zq = my.specs, unif = op.specs,
                  operation = operation,
                  find.cutpoints.zu =
                    .self$fns$find.cutpoints.zu,
                  integralsMult = .self$fns$integralsMult,
                  ExpIntegral = .self$fns$expint)
            },
            stop("ZQuadratic: Invalid operation, stage 2.")
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
            # Normal - ZQuadratic
            tmp.coord <- .self$iget.multiply.constant(-1.0)
            tmp <- ZQuadratic(x = tmp.coord$x, y = tmp.coord$y)
            tmp.specs <- list(
              n = length(tmp$x), x = tmp$x,
              dx = tmp$dx, b0 = tmp$beta0,
              b1 = tmp$beta1, b2 = tmp$beta2,
              normct = tmp$normct)
            coords <-
              .self$fns$convolution.zquadratic.normal(
                zq = tmp.specs, normal = op.specs,
                operation = "+",
                integrals = .self$fns$integrals.quadratic,
                erf = .self$fns$erf, erfi = .self$fns$erfi)
          } else {
            coords <-
              .self$fns$convolution.zquadratic.normal(
                zq = my.specs, normal = op.specs,
                operation = operation,
                integrals = .self$fns$integrals.quadratic,
                erf = .self$fns$erf, erfi = .self$fns$erfi)
          }
        },
        stop("ZQuadratic: Invalid operation, stage 2.")
      )
      return(coords)
    }
  )
)
