#--- Copula tests ------------------------------------------------------

## library(LocalCop)
## library(TMB)
## library(testthat)

library(mvtnorm)

# Generate Values
b1 <- runif(100, -5, 5)
b2 <- runif(100, -5, 5)
rho <- runif(100, 0, 1)

data <- list(model = "pbvn")
parameters <- list(b1 = b1, b2 = b2, rho = rho)

# TMB version of the function
integrated_tmb <- TMB::MakeADFun(data = data, parameters = parameters, silent = TRUE)

pbvn_r <- function(x, b1, b2, rho) {
  loc <- rho * x
  scale <- sqrt(1 - rho * rho)
  ans <- pnorm((b2 - loc) / scale, 0, 1) * dnorm(x, 0, 1)
  ans
}

# R version of function
integrated_r <- function(b1, b2, rho) {
  n <- length(rho)
  ans <- vector(length = n)

  for (i in 1:n) {
    ans[i] <- integrate(function(x) pbvn_r(x, b1[i], b2[i], rho[i]),
                        lower = -Inf, upper = b1[i],
                        subdivisions = 1000L)$value
  }

  return(sum(ans))
}

# mvtnorm version of function
integrated_mvtnorm <- function(b1, b2, rho) {
  n <- length(rho)
  ans <- vector(length = n)

  for (i in 1:n) {
    upper <- c(b1[i], b2[i])
    sigma <- matrix(c(1, rho[i], rho[i], 1), nrow = 2)
    ans[i] <- pmvnorm(lower = c(-Inf, -Inf), upper = upper, sigma = sigma)[1]
  }

  return(sum(ans))
}

# test
test_that("TMB function evaluation is equivalent to analytic R solution", {
  val_r <- integrated_r(b1, b2, rho)
  val_tmb <- integrated_tmb$fn(c(b1, b2, rho))
  expect_equal(val_r, val_tmb)
})

test_that("TMB function evaluation is equivalent to mvtnorm", {
  val_mvtnorm <- integrated_mvtnorm(b1, b2, rho)
  val_tmb <- integrated_tmb$fn(c(b1, b2, rho))
  expect_equal(val_mvtnorm, val_tmb)
})
