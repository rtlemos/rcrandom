test_that("Univariate Normal p.d.f.", {
  mn <- 2
  sd <- 3
  s <- opt(fn = function(x) dnorm(x, mn, sd),
              initial = 0, is.logf = FALSE, verbose = FALSE)
  expect_true(abs(s$mode - mn) < 1e-3)
  expect_true(abs(s$covF[1,1] - sd) < 1e-3)
  expect_true(abs(1 / s$precR[1,1] - sd) < 1e-3)
  expect_true(abs(
    integrate(f = s$pdf, lower = mn - 4 * sd,
              upper = mn + 4 * sd)$value - 1) < 1e-3)
  expect_true(abs(s$cdf(mn + 4 * sd) - 1) < 1e-3)
  expect_true(abs(s$invcdf(0.5) - mn) < 1e-3)
})

test_that("Univariate Gamma p.d.f.", {
  sh <- 2
  rt <- 1/3
  s <- opt(fn = function(x) dgamma(x, sh, rt),
           initial = 1, is.logf = FALSE, verbose = FALSE,
           lower = 0)
  expect_true(abs(s$mode - (sh - 1) / rt) < 1e-3)
  ub <- sh / rt + 6 * sqrt(sh) / rt
  expect_true(abs(
    integrate(f = s$pdf,
              lower = 0,
              upper = ub)$value - 1) < 1e-2)
  expect_true(abs(s$cdf(ub) - 1) < 1e-3)
  expect_true(abs(s$invcdf(1e-50)) < 1e-1)
})

test_that("Univariate Beta p.d.f.", {
  a <- 1.2
  b <- 1.4
  s <- opt(fn = function(x) dbeta(x, a, b),
           initial = 0.1, is.logf = FALSE, verbose = FALSE,
           lower = 0, upper = 1)
  expect_true(abs(s$mode - (a - 1) / (a + b - 2)) < 1e-3) #only true if a,b > 1
  expect_true(abs(integrate(f = s$pdf, lower = 0, upper = 1)$value - 1) < 1e-3)
  expect_true(abs(s$cdf(1) - 1) < 1e-3)
  expect_true(abs(s$invcdf(1e-50)) < 1e-1)
})

test_that("Beta coefficients", {
  x <- c(0.1, 0.15, 0.3, 0.4, 0.55, 0.8, 0.9, 0.97, 0.99999)
  nx <- length(x)
  y <- dbeta(x, 1.6, 0.5, log = TRUE)
  b0 <- b1 <- b2 <- rep(NA, 9 - 1)
  b0[1] <- y[1]
  d <- (x[3] - x[1]) / (x[2] - x[1])
  b1[1] <- d / (d - 1) * (y[2] - y[3] / d ^ 2 - y[1] * (1 - 1 / d ^ 2))
  b2[1] <- y[2] - y[1] - b1[1]
  for (i in 2:(9 - 1)) {
    b0[i] <- y[i]
    b1[i] <- (b1[i - 1] + 2 * b2[i - 1]) * (x[i+1]-x[i]) / (x[i]-x[i-1])
    b2[i] <- y[i + 1] - y[i] - b1[i]
  }
  
  d <- (x[3] - x[1]) / (x[2] - x[1])
  b0 <- y
  b1 <- rep(NA, nx - 1)
  b1[1] <- d / (d - 1) * (y[2] - y[3] / d ^ 2 - y[1] * (1 - 1 / d ^ 2))
  for (i in 2:(nx - 1)) {
    b1[i] <- (b1[i - 1] + 2 * (y[i] - y[i - 1] - b1[i - 1])) * 
      (x[i + 1] - x[i]) / (x[i] - x[i - 1])
  }
  b2[1:(nx - 1)] <- y[2:nx] - y[1:(nx - 1)] - b1[1:(nx - 1)]  
  
  xx <- as.numeric(mapply(1:8, FUN = function(i) {
    seq(x[i], x[i + 1], length = 100)
  }))
  yy <- as.numeric(mapply(1:8, FUN = function(i) {
    myx <- seq(0, 1, length = 100)
    myy <- b0[i] + b1[i] * myx + b2[i] * myx ^ 2
    return(myy)
  }))
  plot(x, y)
  lines(xx, yy)
})
