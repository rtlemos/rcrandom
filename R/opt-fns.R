#' Functions for continuous variables, written in R
#'
#' @return List of R functions for continuous variables
#'
get.rfns <- function(){
  r.functions <- list(
    
    erf = function(x){
      "approximate error function for real values"
      
      t  <-  1.0 / (1.0 + 0.327591100 * abs(x))
      r  <-  1.0 - (0.254829592 * t - 0.284496736 * t ^ 2 +
                      1.421413741 * t ^ 3 -
                      1.453152027 * t ^ 4 +
                      1.061405429 * t ^ 5) * exp(-x ^ 2)
      return(sign(x) * r)
    },
    
    erfi = function(xx){
      "approximate error function for complex values"
      
      u  <- complex(imaginary = abs(xx))
      t  <-  1.0 / (1.0 + 0.327591100 * u )
      r  <-  1.0 - (0.254829592 * t - 0.284496736 * t ^ 2 +
                      1.421413741 * t ^ 3 -
                      1.453152027 * t ^ 4 +
                      1.061405429 * t ^ 5) * exp(-u ^ 2)
      return(sign(xx) * Im(r))
    },
    
    inv.erf = function(x){
      "approximate inverse error function for real values"
      t   <- log(1 - x ^ 2)
      a   <- 0.145
      tpa <- 2.0 / (pi * a)
      r   <- sqrt( sqrt( (tpa + 0.5 * t) ^ 2 - t / a) -
                     (tpa + 0.5 * t))
      return(sign(x) * r)
    },
    
    inv.erfi = function(xx){
      "approximate inverse error function for
      complex values"
      
      u  <- complex(imaginary = abs(xx))
      t  <- log(1.0 - u ^ 2)
      a   <- 0.1382
      tpa <- 2.0 / (pi * a)
      r   <- sqrt( sqrt( (tpa + 0.5 * t) ^ 2 - t / a) -
                     (tpa + 0.5 * t))
      return(sign(xx) * Im(r))
    },
    
    betas.linear = function(nx, x, y){
      "decomposes a vector y into intercept,
      1st and 2nd degree coeffs"
      
      b0 <- y[1:(nx - 1)]
      b1 <- y[2:nx] - y[1:(nx - 1)]
      b2 <- rep(0, nx - 1)
      return(list(b0 = b0, b1 = b1, b2 = b2))
    },
    
    betas.quadratic = function(nx, x, y){
      "same as betas.linear but for ZQuadratic objects"
      
      m <- (nx + 1) / 2 #position of the mode
      b0 <- b1 <- b2 <- rep(NA, nx - 1)
      b0[m] <- y[m]
      b1[m] <- 0.0
      b2[m] <- y[m + 1] - y[m]
      for (i in seq(m + 1, nx - 1)) {
        b0[i] <- y[i]
        b1[i] <- (b1[i - 1] + 2.0 * b2[i - 1]) *
          (x[i + 1] - x[i]) / (x[i] - x[i - 1])
        b2[i] <- y[i + 1] - y[i] - b1[i]
      }
      for (i in seq(m - 1, 1, -1)) {
        b0[i] <- y[i]
        b1[i] <- 2 * (y[i + 1] - y[i]) - b1[i + 1] *
          (x[i + 1] - x[i]) / (x[i + 2] - x[i + 1])
        b2[i] <- y[i + 1] - y[i] - b1[i]
      }
      return(list(b0 = b0, b1 = b1, b2 = b2))
    },
    
    betas.quadratic.old = function(nx, x, y){
      "same as betas.linear but for ZQuadratic objects; old version"
      
      d <- (x[3] - x[1]) / (x[2] - x[1])
      b0 <- y[1:(nx - 1)]
      b1 <- rep(NA, nx - 1)
      b1[1] <- d / (d - 1) * (y[2] - y[3] / d ^ 2 - y[1] * (1 - 1 / d ^ 2))
      for (i in 2:(nx - 1)) {
        b1[i] <- (b1[i - 1] + 2 * (y[i] - y[i - 1] - b1[i - 1])) * 
          (x[i + 1] - x[i]) / (x[i] - x[i - 1])
      }
      b2 <- y[2:nx] - y[1:(nx - 1)] - b1[1:(nx - 1)]      
      return(list(b0 = b0, b1 = b1, b2 = b2))
    },
    
    integrals.linear = function(beta0, beta1, xfrom = 0,
                                xto = 1, ...){
      "integrals for ZLinear objects, used in cdf computation"
      
      out <- exp(beta0) * (exp(beta1 * xto) - exp(beta1 * xfrom)) / beta1
      return(out)
    },
    
    integrals.quadratic = function(beta0, beta1, beta2, xfrom = 0, xto = 1,
                                   erf, erfi){
      "integrals for ZQuadratic objects, used in cdf computation"
      
      twosqpi <- 2 / sqrt(pi)
      spi2    <- sqrt(pi) / 2
      out <- mapply(
        beta0, beta1, beta2, xfrom, xto,
        FUN = function(b0, b1, b2, xfrom, xto){
          imL <- -sign(b2) * (b2 * xfrom + 0.5 * b1) / sqrt(abs(b2))
          imU <- -sign(b2) * (b2 * xto + 0.5 * b1) / sqrt(abs(b2))
          if (b2 == 0) {
            res <- 0
          } else {
            if (abs(b1 / b2) < 100) {
              # do not use approximation
              if (b2 < 0) {
                dlt <- erf(imU) - erf(imL)
              } else {
                dlt <- erfi(imU) - erfi(imL)
              }
            } else if (b2 < 0) {
              # use Taylor series approx. to
              # erf difference (2 terms)
              dlt <- exp(-imL ^ 2) * twosqpi * (imU - imL) *
                (1 - (imU - imL) * imL)
            } else {
              # use T.S. approximation to
              # erfi difference (2 terms)
              dlt <- exp(imL ^ 2) * twosqpi * (imU - imL) *
                (1 + (imU - imL) * imL)
            }
            res <- -spi2 / sqrt(abs(b2)) * dlt *
              exp((4 * b0 * b2 - b1 ^ 2) / (4 * b2))
          }
          return(res)
        })
      out[beta2 == 0] <- 0
      return(out)
    },
    
    integralsMult = function(beta0, beta1, beta2, xfrom, xto, ExpIntegral){
      "to be documented"
      
      out <- mapply(beta0, beta1, beta2, xfrom, xto,
        FUN = function(beta0, beta1, beta2, xfrom, xto){
          if (sign(xfrom) == sign(xto)) {
            xmiddle1 <- xmiddle2 <- 0.5 * (xfrom + xto)
          } else {
            xmiddle1 <- -1e-16
            xmiddle2 <- +1e-16
          }
          xL <- c(xfrom, xmiddle1)
          xU <- c(xmiddle2, xto)
          integr <- mapply(1:2, FUN = function(i){
            kappa0 <- beta0 - beta2 * (xL[i] * xU[i])
            kappa1 <- beta1 + beta2 * (xL[i] + xU[i])
            eL <- sign(xL[i]) *
              ExpIntegral(kappa1 * xL[i])
            eU <- sign(xU[i]) *
              ExpIntegral(kappa1 * xU[i])
            integral <- exp(kappa0) * (eU - eL)
            return(integral)
          })
          res <- sum(integr)
          return(res)
        })
      return(out)
    },
    
    pdf.linear = function(z, x, dx, beta0, beta1, normct,
                          pdf.quadratic, do.log){
      "compute the probability density function for
      a ZLinear object"
      
      beta2 <- rep(0,length(beta1))
      return(pdf.quadratic(z, x, dx, beta0, beta1, beta2,
                           normct, do.log))
    },
    
    pdf.quadratic = function(z, x, dx, beta0, beta1,
                             beta2, normct, do.log){
      "compute the probability density function for
      a ZQuadratic object"
      
      pos <- mapply(z, FUN = function(z){
        j   <- which.min(abs(z - x))
        out <- if (x[j] > z & j > 1) j - 1 else j
        return(out)
      })
      zscaled <- (z - x[pos]) / dx[pos]
      logpdf <- beta0[pos] + beta1[pos] * zscaled +
        beta2[pos] * zscaled ^ 2
      if (do.log) {
        out <- logpdf - log(normct)
      } else {
        out <- exp(logpdf) / normct
      }
      out[is.na(out)] <- if (do.log) -50 else 0.0
      return(out)
    },
    
    cdf.linear = function(z, x, dx, beta0, beta1, pcdf, normct){
      "compute the cumulative density function for
      a ZLinear object"
      
      pos <- mapply(z, FUN = function(z){
        j   <- which.min(abs(z - x))
        out <- if (x[j] > z) j - 1 else j
        return(out)
      })
      out <- rep(0,length(z))
      n   <- length(x)
      for (i in 1:(n - 1)) {
        b0 <- beta0[i]
        b1 <- beta1[i]
        zp <- (z[pos == i] - x[i]) / dx[i]
        if (length(zp) > 0) {
          cp  <- exp(b0) * (exp(b1 * zp) - 1.0) / b1
          rd  <- if (i > 1) pcdf[i - 1] else 0.0
          r 	<- rd + cp * dx[i] / normct
          out[pos == i] <- r
        }
      }
      out[pos >= n] <- 1
      return(out)
    },
    
    invcdf.linear = function(p, x, dx, beta0, beta1, pcdf, normct){
      "compute the inverse CDF for a ZLinear object"
      
      pos <- mapply(p, FUN = function(p){
        j <- which.min(abs(p - pcdf))
        posi <- if (pcdf[j] > p & j > 1) j - 1 else j
        return(posi)
      })
      out <- rep(0,length(p))
      n   <- length(x)
      for (i in 1:(n - 1)) if (length(p[pos == i]) > 0) {
        b0 <- beta0[i]
        b1 <- beta1[i]
        rd <- if (i > 1) pcdf[i - 1] else 0.0
        zp <- log((p[pos == i] - rd) * normct * b1 /
                    (dx[i] * exp(b0)) + 1.0) / b1
        r 	<- x[i] + zp * dx[i]
        out[pos == i] <- r
      }
      return(out)
    },
    
    cdf.quadratic = function(z, x, dx, beta0, beta1, beta2, pcdf,
                             normct, erf, erfi){
      "compute the cumulative density function for a ZQuadratic object"
      
      pos <- mapply(z, FUN = function(z){
        j <- which.min(abs(z - x))
        out <- if (x[j] > z & j > 1) j - 1 else j
        return(out)
      })
      out <- rep(0,length(z))
      twosqpi <- 2/sqrt(pi)
      spi2 <- sqrt(pi)/2
      n <- length(x)
      for (i in 1:(n - 1)) if (sum(pos == i) > 0) {
        cr <- (pos == i)
        b0 <- beta0[i]
        b1 <- beta1[i]
        b2 <- beta2[i]
        zp <- (z[cr] - x[i]) / dx[i]
        sb2 <- sqrt(abs(b2))
        imL <- -sign(b2) * 0.5 * b1 / sb2
        imU <- -sign(b2) * (b2 * zp + 0.5 * b1) / sb2
        if (b2 == 0) {
          r <- 0
        } else {
          if (abs(b1 / b2) < 100) { #do not use approximation
            if (b2 < 0) {
              dlt <- erf(imU) - erf(imL)
            } else {
              dlt <- erfi(imU) - erfi(imL)
            }
          } else if (b2 < 0) {
            # 2-term Taylor series approx. to erf difference
            dlt <- exp(-imL ^ 2) * twosqpi * (imU - imL) *
              (1 - (imU - imL) * imL)
          } else {
            # 2-term T.S. approximation to erfi difference
            dlt <- exp(imL ^ 2) * twosqpi * (imU - imL) *
              (1 + (imU - imL) * imL)
          }
          cp <- -spi2 / sb2 * dlt *
            exp((4 * b0 * b2 - b1 ^ 2) / (4 * b2))
          rd <- if (i > 1) pcdf[i - 1] else 0.0
          r  <- rd + cp * dx[i] / normct
        }
        out[cr] <- r
      }
      out[pos >= n] <- 1
      return(out)
    },
    
    invcdf.quadratic = function(p, x, dx, beta0, beta1,
                                beta2, pcdf, normct,
                                erf, erfi, inv.erf,
                                inv.erfi){
      "compute the inverse CDF for a ZQuadratic object"
      
      n   <- length(x)
      pos <- mapply(p, FUN = function(p){
        j    <- which.min(abs(p - pcdf))
        posi <- if (pcdf[j] < p & j < n) j + 1 else j
        return(posi)
      })
      out     <- rep(0, length(p))
      twosqpi <- 2 / sqrt(pi)
      spi2    <- sqrt(pi) / 2
      for (i in 1:(n - 1)) if (sum(pos == i) > 0) {
        cr <- (pos == i)
        b0  <- beta0[i]
        b1  <- beta1[i]
        b2  <- beta2[i]
        xi  <- x[i]
        dxi <- dx[i]
        rd  <- if (i > 1) pcdf[i - 1] else 0.0
        sb2 <- sqrt(abs(b2))
        cp  <- (normct / dxi) * (p[cr] - rd)
        dlt <- -twosqpi * cp * sb2 /
          exp( (4 * b0 * b2 - b1 ^ 2) / (4 * b2))
        imL <- -sign(b2) * 0.5 * b1 / sb2
        if ( abs(b1 / b2) < 100) {
          #do not use approximation
          if (b2 < 0) {
            imU <- inv.erf(dlt + erf(imL))
          } else {
            imU <- inv.erfi(dlt + erfi(imL))
          }
          zp  <- -imU / sb2 - 0.5 * b1 / b2
        } else {
          if (b2 < 0) {
            # 2-term Taylor series approx. to erf difference
            a <- imL
            b <- -1 - 2 * imL ^ 2
            c <- imL ^ 3 + imL + dlt * spi2 / exp(-imL ^ 2)
          } else {
            # 2-term T.S. approximation to erfi difference
            a <- imL
            b <- 1 - 2 * imL ^ 2
            c <- imL ^ 3 - imL - dlt * spi2 / exp(imL ^ 2)
          }
          root1 <- (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
          root2 <- (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
          zp  <- mapply(root1, root2,
                        FUN = function(r1, r2){
                          z1 <- -r1 / sb2 - 0.5 * b1 / b2
                          z2 <- -r2 / sb2 - 0.5 * b1 / b2
                          if (z1 >= 0 & z1 <= 1) {
                            res <- z1
                          } else if (z2 >= 0 & z2 <= 1) {
                            res <- z2
                          } else {
                            stop("zp should be between 0 and 1!")
                          }
                          return(res)
                        })
        }
        out[cr] <- xi + dxi * zp
      }
      out[pos >= n] <- x[n]
      return(out)
    },
    
    convolution.zquadratic.zquadratic = function(
      l, r, operation, find.cutpoints.zz,
      integrals, erf, erfi){
      "performs a convolution of two ZQuadratic objects"
      
      if (operation != "+") stop("not implemented")
      nl <- l$n
      xl <- l$x
      dxl <- l$dx
      b0l <- l$b0
      b1l <- l$b1
      b2l <- l$b2
      normctl <- l$normct
      nr <- r$n
      xr <- r$x
      dxr <- r$dx
      b0r <- r$b0
      b1r <- r$b1
      b2r <- r$b2
      normctr <- r$normct
      xlol <- xl[1:(nl - 1)]
      xhil <- xl[2:nl]
      xlor <- xr[1:(nr - 1)]
      xhir <- xr[2:nr]
      a0l <- b0l - log(normctl)
      a0r <- b0r - log(normctr) - b1r * xlor / dxr +
        b2r * (xlor / dxr) ^ 2
      a1l <- -b1l/dxl
      a1r <- b1r/dxr - 2 * b2r * xlor / dxr ^ 2
      a2l <- b2l / dxl ^ 2
      a2r <- b2r / dxr ^ 2
      st <- xl[1]  + xr[1] + 1e-9
      en <- xl[nl] + xr[nr] - 1e-9
      n <- 0.5*(nr + nl)
      x <- y <- rep(NA, n)
      for (i in 1:n) {
        evaluation.point <- st + (en - st) *
          (i - 1) / (n - 1)
        cutp <- find.cutpoints.zz(nl, xl, xlol, xhil, nr,
                                  xr, xlor, xhir,
                                  s = evaluation.point)
        if (!is.null(cutp)) {
          ylo  <- cutp$ylo
          yhi  <- cutp$yhi
          il   <- cutp$indl
          ir   <- cutp$indr
          alpha0 <- a0r[ir] + a0l[il] + b1l[il] *
            (evaluation.point - xlol[il]) / dxl[il] +
            b2l[il] *
            ((evaluation.point - xlol[il]) / dxl[il]) ^ 2
          alpha1 <- a1r[ir] + a1l[il] - 2 * b2l[il] *
            (evaluation.point - xlol[il]) / dxl[il] ^ 2
          alpha2 <- a2r[ir] + a2l[il]
          ints <- integrals(alpha0, alpha1, alpha2,
                            xfrom = ylo, xto = yhi,
                            erf = erf, erfi = erfi)
          x[i] <- evaluation.point
          y[i] <- log(sum(ints))
        } else {
          stop("Empty vector of cutpoints.")
        }
      }
      return(list(x = x, y = y))
    },
    
    find.cutpoints.zz = function(nl = NULL, xl, xlol = NULL,
                                 xhil = NULL, nr = NULL, xr,
                                 xlor = NULL, xhir = NULL,
                                 s){
      "to be documented"
      
      nl <- if (is.null(nl)) length(xl) else nl
      nr <- if (is.null(nr)) length(xr) else nr
      xlol <- if (is.null(xlol)) xl[1:(nl - 1)] else xlol
      xlor <- if (is.null(xlor)) xr[1:(nr - 1)] else xlor
      xhil <- if (is.null(xhil)) xl[2:nl] else xhil
      xhir <- if (is.null(xhir)) xr[2:nr] else xhir
      xml  <- 0.5 * (xl[1:(nl - 1)] + xl[2:nl] )
      xmr  <- 0.5 * (xr[1:(nr - 1)] + xr[2:nr] )
      
      ylo  <- yhi <- indl <- indr <- rep(NA, nl + nr)
      lo   <- max(s - xl[nl], xr[1])
      hiF  <- min(s - xl[ 1], xr[nr])
      i    <- 0
      while (lo < hiF) {
        i    <- i + 1
        locl <- which(xlol < (s - lo) & (s - lo) <= xhil)
        if (locl > 1 & abs(s - lo - xlol[locl]) < 1e-10 ) {
          locl <- locl - 1
        }
        locr <- which(xlor <= lo & lo <  xhir)
        if (locr < nr & abs( xhir[locr] - lo) < 1e-10 ) {
          locr <- locr + 1
        }
        dl   <- (s - lo) - xlol[locl]
        dr   <- xhir[locr] - lo
        df   <- hiF - lo
        hi   <- lo + min(dl, dr, df)
        ylo[i]  <- lo
        yhi[i]  <- hi
        indl[i] <- locl
        indr[i] <- locr
        if (i > nl + nr) {
          stop("error: algorithm has stalled")
        }
        lo	 <- hi
      }
      if (i > 0) {
        return(list(ylo = ylo[1:i], yhi = yhi[1:i],
                    indl = indl[1:i], indr = indr[1:i]))
      } else {
        return(NULL)
      }
    },
    
    find.cutpoints.zu = function(n = NULL, x, xlow = NULL,
                                 xhi = NULL, unif.lb,
                                 unif.ub, s, operation){
      "to be documented"
      
      n    <- if (is.null(n)) length(x) else n
      xlow <- if (is.null(xlow)) x[1:(n - 1)] else xlow
      xhi  <- if (is.null(xhi))  x[2:n] else xhi
      xm   <- 0.5 * (x[1:(n - 1)] + x[2:n])
      ylo  <- yhi <- ind <- rep(NA, n)
      switch(operation,
             "+" = {aux  <- s - c(unif.ub, unif.lb)},
             "-" = {aux  <- s + c(unif.lb, unif.ub)},
             "*" = {rati <- s / c(unif.ub, unif.lb)
             if (rati[1] < rati[2]) {
               aux  <- rati
             } else {
               aux <- rev(rati)
             }},
             stop("find.cutpoints.zu: invalid operation")
      )
      lo <- max(x[1], aux[1])
      hiF <- min(x[n], aux[2])
      step <- (hiF - lo) / (n - 1)
      if (abs(lo - hiF) <= 1e-10) {
        i <- 1
        ind[i] <- which( xlow <=  lo & lo < xhi)
        ylo[i] <- lo
        yhi[i] <- hiF
      } else {
        i <- 0
        while (abs(lo - hiF) > 1e-10) {
          i <- i + 1
          loc <- which( xlow <=  lo & lo < xhi)
          if (loc > 1  & abs(lo - xlow[loc]) < 1e-10 ) {
            loc <- loc - 1
          }
          #dl   <- w - xlow[locl]
          #df   <- hiF - lo
          hi     <- lo + step
          ylo[i] <- lo
          yhi[i] <- hi
          ind[i] <- loc
          if (i > n) {
            print(c(lo, hiF))
            stop("error: algorithm has stalled")
          }
          lo	 <- hi
        }
      }
      if (i > 0) {
        out <- list(ylo = ylo[1:i], yhi = yhi[1:i],
                    ind = ind[1:i])
      } else {
        out <- NULL
      }
      return(out)
    },
    
    convolution.zquadratic.normal = function(
      zq, normal, operation, integrals, erf = NULL,
      erfi = NULL){
      "performs a convolution of ZQuadratic and
      Normal objects"
      
      nl <- zq$n
      xl <- zq$x
      dxl <- zq$dx
      b0l <- zq$b0
      b1l <- zq$b1
      b2l <- zq$b2
      normctl <- zq$normct
      m  <- normal$m; v <- normal$v
      n <- length(xl)
      x <- y <- rep(NA, n)
      ints   <- matrix(nrow = n - 1, ncol = n)
      const  <- -0.5 * log(2 * pi * v) - log(normctl)
      xlow   <- xl[1:(n - 1)]
      xhi    <- xl[2:n]
      a0part <- b0l - b1l * xlow / dxl +
        b2l * (xlow / dxl) ^ 2 + const
      a1part <- b1l/dxl - 2 * b2l * xlow / dxl ^ 2
      alpha2 <- -0.5 / v + b2l / dxl ^ 2
      for (i in 1:n) {
        evaluation.point <- get(operation)(xl[i], m) -
          3.5 * sqrt(v) + 7.0 * (i - 1) * sqrt(v) / (n - 1)
        alpha0 <- a0part -
          0.5 * (evaluation.point - m) ^ 2 / v
        alpha1 <- get(operation)(a1part,
                                 (evaluation.point - m) / v)
        ints[,i] <- integrals(alpha0, alpha1, alpha2,
                              xfrom = xlow, xto = xhi,
                              erf = erf, erfi = erfi)
        x[i] <- evaluation.point
        y[i] <- log(sum(ints[, i]))
      }
      return(list(x = x, y = y))
    },
    
    convolution.zquadratic.uniform = function(
      zq, unif, operation, find.cutpoints.zu,
      integralsMult, ExpIntegral){
      "performs a convolution of ZQuadratic and
      Uniform objects"
      
      nl <- zq$n
      xl <- zq$x
      dxl <- zq$dx
      b0l <- zq$b0
      b1l <- zq$b1
      b2l <- zq$b2
      normctl <- zq$normct
      unif.lb <- unif$lb
      unif.ub <- unif$ub
      unif.n <- unif$n
      n <- length(xl)
      x <- y <- rep(NA, n)
      const  <- -log(normctl * (unif.ub - unif.lb))
      xlow   <- xl[1:(n - 1)]
      xhi    <- xl[2:n]
      a0     <- b0l - b1l * xlow / dxl +
        b2l * (xlow / dxl) ^ 2 + const
      a1     <- b1l / dxl - 2 * b2l * xlow / dxl ^ 2
      a2     <- b2l / dxl ^ 2
      switch(
        operation,
        "+" = {cn <- c(xl[1] + unif.lb, xl[n] + unif.ub)},
        "-" = {cn <- c(xl[1] - unif.ub, xl[n] - unif.lb)},
        "*" = {cn <- c(xl[1] * unif.lb, xl[1] * unif.ub,
                       xl[n] * unif.lb, xl[n] * unif.ub)},
        stop("Invalid operation.")
      )
      st <- min(cn) + 1e-14
      rn <- max(cn) - min(cn) - 2e-14
      for (i in 1:n) {
        candidate <- st + rn * (i - 1) / (n - 1)
        if (abs(candidate) > 1e-16) {
          evaluation.point <- candidate
        } else {
          evaluation.point <- 1e-16
        }
        cutp <- find.cutpoints.zu(
          n, xl, xlow, xhi, unif.lb,
          unif.ub, s = evaluation.point,
          operation = operation)
        if (!is.null(cutp)) {
          ylo    <- cutp$ylo
          yhi    <- cutp$yhi
          il     <- cutp$ind
          alpha0 <- a0[il]
          alpha1 <- a1[il]
          alpha2 <- a2[il]
          ints <- integralsMult(alpha0, alpha1, alpha2,
                                xfrom = ylo, xto = yhi,
                                ExpIntegral = ExpIntegral)
          x[i] <- evaluation.point
          y[i] <- log(sum(ints))
        } else {
          print(cutp)
          stop("Empty vector of cutpoints")
        }
      }
      return(list(x = x, y = y))
    },
    
    expint = function(x){
      "The exponential integral function"
      
      if (x == 0) {
        return(-Inf)
      } else {
        out <- .Fortran("expint", x = as.double(x),
                        res = as.double(0))
        return(out$res)
      }
    }
  )
  return(r.functions)
  }

#' Functions for continuous variables, written in Fortran
#'
#' @return List of functions Functions for continuous
#' variables, in Fortran
#'
get.ffns <- function(){
  fortran.functions <- list(
    
    erf = function(x){
      "approximate error function for real values"
      
      nx  <- length(x)
      out <- .Fortran("erfV", n = as.integer(nx),
                      x = as.double(x),
                      res = as.double(rep(0, nx)))
      return(out$res)
    },
    
    erfi = function(x){
      "approximate error function for complex values"
      
      nx  <- length(x)
      out <- .Fortran("erfiV", n = as.integer(nx),
                      x = as.double(x),
                      res = as.double(rep(0, nx)))
      return(out$res)
    },
    
    inv.erf = function(x){
      "approximate inverse error function for real values"
      
      nx  <- length(x)
      out <- .Fortran("inv_erfV", n = as.integer(nx),
                      x = as.double(x),
                      res = as.double(rep(0, nx)))
      return(out$res)
    },
    
    inv.erfi = function(x){
      "approximate inverse error function for
      complex values"
      
      nx  <- length(x)
      out <- .Fortran("inv_erfiV", n = as.integer(nx),
                      x = as.double(x),
                      res = as.double(rep(0, nx)))
      return(out$res)
    },
    
    betas.linear = function(nx, x, y){
      "decomposes a vector y into intercept,
      1st and 2nd degree coeffs"
      
      emp <- rep(0, nx - 1)
      out <- .Fortran("continuous_betas_linear",
                      npt = as.integer(nx),
                      y = as.double(y),
                      b0 = as.double(emp),
                      b1 = as.double(emp))
      return(list(b0 = out$b0, b1 = out$b1, b2 = emp))
    },
    
    betas.quadratic = function(nx, x, y){
      "same as betas.linear but for ZQuadratic objects"
      
      emp <- rep(0, nx - 1)
      out <- .Fortran("continuous_betas_quadratic",
                      npt = as.integer(nx),
                      x = as.double(x),
                      y = as.double(y),
                      b0 = as.double(emp),
                      b1 = as.double(emp),
                      b2 = as.double(emp))
      return(list(b0 = out$b0, b1 = out$b1, b2 = out$b2))
    },
    
    integrals.linear = function(beta0, beta1, xfrom = 0,
                                xto = 1){
      "integrals for ZLinear objects,
      used in cdf computation"
      
      if (length(xfrom) == 1) xfrom <- rep(xfrom,
                                           length(beta0))
      if (length(xto)   == 1) xto   <- rep(xto,
                                           length(beta0))
      out <- .Fortran("continuous_integrals_b2equal0",
                      npt = as.integer(length(beta0) + 1),
                      beta0 = as.double(beta0),
                      beta1 = as.double(beta1),
                      xlow = as.double(xfrom),
                      xhigh = as.double(xto),
                      integrals = as.double(
                        rep(0, length(beta0))))
      return(out$integrals)
    },
    
    integrals.quadratic = function(beta0, beta1, beta2,
                                   xlow, xhigh, ...){
      "integrals for ZQuadratic objects, used in cdf computation"
      
      if (length(xlow) == 1) xlow = rep(xlow, length(beta0))
      if (length(xhigh) == 1) xhigh = rep(xhigh,
                                          length(beta0))
      out <- .Fortran("continuous_integrals",
                      npt = as.integer(length(beta0) + 1),
                      beta0 = as.double(beta0),
                      beta1 = as.double(beta1),
                      beta2 = as.double(beta2),
                      xlow = as.double(xlow),
                      xhigh = as.double(xhigh),
                      integrals = as.double(
                        rep(0, length(beta0))))
      return(out$integrals)
    },
    
    pdf.linear = function(z, x, dx, beta0, beta1, normct,
                          do.log, pdf.quadratic){
      "compute the probability density function for
      a ZLinear object"
      
      beta2 <- rep(0, length(beta1))
      out <- pdf.quadratic(z, x, dx, beta0, beta1, beta2,
                           do.log, normct)
      return(out)
    },
    
    pdf.quadratic = function(z, x, dx, beta0, beta1, beta2,
                             normct, do.log){
      "compute the probability density function for a
      ZQuadratic object"
      
      out <- .Fortran("continuous_pdf",
                      n = as.integer(length(z)),
                      z = as.double(z),
                      npt = as.integer(length(x)),
                      x = as.double(x), dx = as.double(dx),
                      beta0 = as.double(beta0),
                      beta1 = as.double(beta1),
                      beta2 = as.double(beta2),
                      normct = as.double(normct),
                      dolog = do.log,
                      pdf = as.double(rep(0, length(z))))
      return(out$pdf)
    },
    
    cdf.linear = function(z, x, dx, beta0, beta1, pcdf,
                          normct, ...){
      "compute the cumulative density function for a
      ZLinear object"
      
      out <- .Fortran("continuous_cdf_linear",
                      n = as.integer(length(z)),
                      z = as.double(z),
                      npt = as.integer(length(x)),
                      x = as.double(x), dx = as.double(dx),
                      beta0 = as.double(beta0),
                      beta1 = as.double(beta1),
                      pcdf = as.double(pcdf),
                      normct = as.double(normct),
                      p = as.double(rep(0,length(z))))
      return(out$p)
    },
    
    cdf.quadratic = function(z, x, dx, beta0, beta1, beta2,
                             pcdf, normct, ...){
      "compute the cumulative density function for a
      ZQuadratic object"
      
      out <- .Fortran("continuous_cdf_quadratic",
                      n = as.integer(length(z)),
                      z = as.double(z),
                      npt = as.integer(length(x)),
                      x = as.double(x), dx = as.double(dx),
                      beta0 = as.double(beta0),
                      beta1 = as.double(beta1),
                      beta2 = as.double(beta2),
                      pcdf = as.double(pcdf),
                      normct = as.double(normct),
                      p = as.double(rep(0, length(z))))
      return(out$p)
    },
    
    invcdf.linear = function(p, x, dx, beta0, beta1,
                             pcdf, normct){
      "compute the inverse CDF for a ZLinear object"
      
      out <- .Fortran("continuous_invcdf_linear",
                      n = as.integer(length(p)),
                      p = as.double(p),
                      npt = as.integer(length(x)),
                      x = as.double(x), dx = as.double(dx),
                      beta0 = as.double(beta0),
                      beta1 = as.double(beta1),
                      pcdf = as.double(pcdf),
                      normct = as.double(normct),
                      z = as.double(rep(0, length(p))))
      return(out$z)
    },
    
    invcdf.quadratic = function(p, x, dx, beta0, beta1,
                                beta2, pcdf, normct, erf,
                                erfi, inv.erf, inv.erfi){
      "compute the inverse CDF for a ZQuadratic object"
      
      out <- .Fortran("continuous_invcdf_quadratic",
                      n = as.integer(length(p)),
                      p = as.double(p),
                      npt = as.integer(length(x)),
                      x = as.double(x), dx = as.double(dx),
                      beta0 = as.double(beta0),
                      beta1 = as.double(beta1),
                      beta2 = as.double(beta2),
                      pcdf = as.double(pcdf),
                      normct = as.double(normct),
                      z = as.double(rep(0, length(p))))
      return(out$z)
    },
    
    convolution.zquadratic.zquadratic = function(
      l, r, operation, ...){
      "performs a convolution of two ZQuadratic objects"
      
      if (operation != "+") stop("not implemented")
      nl <- l$n
      xl <- l$x
      dxl <- l$dx
      b0l <- l$b0
      b1l <- l$b1
      b2l <- l$b2
      normctl <- l$normct
      nr <- r$n
      xr <- r$x
      dxr <- r$dx
      b0r <- r$b0
      b1r <- r$b1
      b2r <- r$b2
      normctr <- r$normct
      emp <- rep(0, (nl + nr) / 2)
      stop("To be coded in Fortran incl. find.cutpoints.zz")
      out <- .Fortran("convolution.zquadratic.zquadratic",
                      nl = as.integer(nl),
                      xl = as.double(xl),
                      dxl = as.double(dxl),
                      b0l = as.double(b0l),
                      b1l = as.double(b1l),
                      b2l = as.double(b2l),
                      nr = as.integer(nr),
                      xr = as.double(xr),
                      dxr = as.double(dxr),
                      b0r = as.double(b0r),
                      b1r = as.double(b1r),
                      b2r = as.double(b2r),
                      x = as.double(emp),
                      y = as.double(emp))
      return(list(x = out$x, y = out$y))
    },
    
    convolution.zquadratic.normal = function(
      xl, dxl, b0l, b1l, b2l, normct, m, v){
      "performs a convolution of ZQuadratic and
      Normal objects"
      
      stop("To Be coded in Fortran")
      nl  <- length(xl)
      emp <- rep(0, nl)
      out <- .Fortran("convolution_zquadratic_normal",
                      nl = as.integer(nl),
                      xl = as.double(xl),
                      dxl = as.double(dxl),
                      b0l = as.double(b0l),
                      b1l = as.double(b1l),
                      b2l = as.double(b2l),
                      m = as.double(m),
                      v = as.double(v),
                      x = as.double(emp),
                      y = as.double(emp))
      return(list(x = out$x, y = out$y))
    },
    
    convolution.zquadratic.uniform = function(
      zq, unif, operation, find.cutpoints.zu){
      "performs a convolution of ZQuadratic and
      Uniform objects"
      
      nl <- zq$n
      xl <- zq$x
      dxl <- zq$dx
      b0l <- zq$b0
      b1l <- zq$b1
      b2l <- zq$b2
      normctl <- zq$normct
      unif.lb <- unif$lb
      unif.ub <- unif$ub
      unif.n <- unif$n
      stop("To Be coded in Fortran")
      nl  <- length(xl)
      emp <- rep(0, nl)
      out <- .Fortran("convolution_zquadratic_uniform",
                      nl = as.integer(nl),
                      xl = as.double(xl),
                      dxl = as.double(dxl),
                      b0l = as.double(b0l),
                      b1l = as.double(b1l),
                      b2l = as.double(b2l),
                      unif_lb = as.double(unif.lb),
                      unif_ub = as.double(unif.ub),
                      x = as.double(emp),
                      y = as.double(emp))
      return(list(x = out$x, y = out$y))
    },
    
    expint = function(x){
      "The exponential integral function"
      
      out <- .Fortran("expint", x = as.double(x),
                      res = as.double(0))
      return(out$res)
    }
  )
  return(fortran.functions)
}
