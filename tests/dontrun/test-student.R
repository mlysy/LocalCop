#--- tests for student-t copula ------------------------------------------------

library(mvtnorm)
library(testthat)

#' Simplified form of bivariate student-t log-pdf.
#'
#' @param y1,y2 Quantiles.
#' @param rho Correlation parameter.
#' @param nu Degrees of freedom parameter.
#' @return Log-pdf of bivariate student-t.
dstudent <- function(y1, y2, rho, nu) {
  Sigma_det <- 1 - rho^2
  ans <- (y1*y1 + y2*y2 - 2*rho * y1*y2) / Sigma_det
  ans <- -(nu + 2)/2 * log(1 + ans/nu)
  ans <- ans - 2 * log(sqrt(2*pi))
  ans <- ans - .5 * log(Sigma_det)
  ans
}

#' Unsimplified version using [mvtnorm::dmvt()].
#'
dstudent_uns <- function(y1, y2, rho, nu) {
  Sigma <- rho + diag(rep(1-rho, 2))
  mvtnorm::dmvt(x = c(y1, y2), sigma = Sigma, df = nu, log = TRUE)
}

n_sim <- 10
y1 <- rnorm(n_sim)
y2 <- rnorm(n_sim)
rho <- runif(n_sim, -1, 1)
nu <- runif(n_sim, 0, 50)

ans_simp <- mapply(dstudent, y1 = y1, y2 = y2, rho = rho, nu = nu)
ans_uns <- mapply(dstudent_uns, y1 = y1, y2 = y2, rho = rho, nu = nu)

expect_equal(ans_simp, ans_uns)
