#' ---
#' title: CDF of bivariate normal and student-t distributions
#' author: Martin Lysy
#' ---
#'
#' # Synopsis
#'
#' Implementation of the bivariate normal and student-t CDF using the one-dimensional integral underlying algorithm QRSVT (above equation 2.9) of Genz and Bretz (2002).  This formula is used rather than the integral underlying algorithm QRSVN (above equation 2.11), which the authors report to be faster in high dimensions, because (i) the performance of the two is very similar in two dimensions and (ii) QRSVT is a one-dimensional integral for both the normal and the student-t, whereas QRSVN is two-dimensional for the latter.

library(mvtnorm)
library(testthat)

#' Calculate the integrand in algorithm QRSVT.
#'
#' Integrating this function gives the rectangular probabilities of a standard bivariate normal or student-t distribution, i.e., which are marginally `N(0, 1)` or `t_(nu)`.
#'
#' @param u Vector of value between 0 and 1.  Corresponds to `w_1` in equation 2.9.
#' @param lower Vector of two lower integration limits.  Correspond to `(a_1, a_2)` in equation 2.9.
#' @param upper Vector of two upper integration limits.  Correspond to `(b_1, b_2)` in equation 2.9.
#' @param rho Correlation parameter.
#' @param nu Degrees of freedom parameter.  Setting to zero means normal CDF.
#'
#' @return Vector of same length as `u` corresponding to the integrand of the QRSVT method.
#'
#' @details  Calculates
#' ```
#' Pr(a_1 < X_1 < b_1, a_2 < X_2 < b_2)
#' ```
#' via the SVT integral with the following symmetries:
#' - If `a_1 > 0`, calculates `Pr(-b_1 < -X_1 < -a_1, a_2 < X_2 < b_2)`, which is identical to the above, but has greater precision given that `pnorm()` and `pt()` are more accurate in the lower tails.
#' - Similarly for `a_2 > 0`.
#' - Orders the variables such that the one with the smaller integration length comes first.
svt <- function(u, lower, upper, rho, nu) {
  lims <- rbind(lower, upper)
  # pnorm and pt are more accurate in lower tails.
  # so flip integral if needed
  if(lims[1,1] > 0) {
    lims[,1] <- -lims[2:1,1]
    rho <- -rho
  }
  if(lims[1,2] > 0) {
    lims[,2] <- -lims[2:1,2]
    rho <- -rho
  }
  # smaller integral first
  lims <- lims[,order(lims[2,] - lims[1,])]
  .svt(u, lower = lims[1,], upper = lims[2,],
       rho = rho, nu = nu)
}

#' Raw calculation.
#'
#' Argument formatting is done in `svt()`.
#' @noRd
.svt <- function(u, lower, upper, rho, nu) {
  if(nu == 0) {
    pfun <- function(q, nu) pnorm(q)
    qfun <- function(p, nu) qnorm(p)
    tfun <- function(y, nu) 1
  } else {
    pfun <- function(q, nu) pt(q, df = nu)
    qfun <- function(p, nu) qt(p, df = nu)
    tfun <- function(y, nu) sqrt((nu+1)/(nu+y^2))
  }
  # a/b is a_tilde/b_tilde in the notation of 2.9
  a1 <- lower[1]
  b1 <- upper[1]
  c12 <- rho
  c22 <- sqrt(1-rho^2)
  d1 <- pfun(a1, nu)
  e1 <- pfun(b1, nu)
  if(e1 - d1 < 1e-15) return(e1 - d1)
  z1 <- d1 + (e1 - d1) * u # exp(log(d1) + log(e1 / d1 - 1) + u
  y1 <- qfun(z1, nu)
  t2 <- tfun(y1, nu)
  a2 <- t2 * (lower[2] - c12 * y1)/c22
  b2 <- t2 * (upper[2] - c12 * y1)/c22
  d2 <- pfun(a2, nu+1)
  e2 <- pfun(b2, nu+1)
  (e1 - d1) * (e2 - d2)
}

#--- tests ---------------------------------------------------------------------

ilogit <- function(x) 1/(1 + exp(-x))
logit <- function(x) log(x) - log(1-x)

test_matrix <- expand.grid(
  has_df = c(FALSE, TRUE),
  lower_inf = c(FALSE, TRUE),
  upper_inf = c(FALSE, TRUE),
  extreme_lim = c(FALSE, TRUE),
  extreme_corr = c(FALSE, TRUE)
)
n_test <- nrow(test_matrix)
n_test_per_case <- 50

for(i in 1:n_test) {
  has_df <- test_matrix$has_df[i]
  lower_inf <- test_matrix$lower_inf[i]
  upper_inf <- test_matrix$upper_inf[i]
  extreme_lim <- test_matrix$extreme_lim[i]
  extreme_corr <- test_matrix$extreme_corr[i]
  sigma <- if(extreme_lim) 20 else 3
  if(has_df) {
    qfun <- function(p, nu) qt(p, nu)
  } else {
    qfun <- function(p, nu) qnorm(p)
  }
  for(j in 1:n_test_per_case) {
    ## u <- sort(runif(1e4))
    rho <- if(extreme_corr) runif(1, .99, 1) else runif(1)
    rho <- sample(c(-1, 1), 1) * rho
    nu <- if(has_df) sample(10, 1) else 0
    lim1 <- sort(rnorm(2, sd = sigma))
    lim2 <- sort(rnorm(2, sd = sigma))
    lower <- c(lim1[1], lim2[1])
    upper <- c(lim1[2], lim2[2])
    if(lower_inf) lower[sample(2, 1)] <- -Inf
    lower <- c(-Inf, -Inf)
    if(upper_inf) upper[sample(2, 1)] <- Inf
    corr <- diag(rep(1-rho, 2)) + rho
    ## p1 <- integrate(function(x) svt(x, lower, upper, rho, nu),
    ##                 lower = 0, upper = 1, subdivisions = 1000L)
    ## p1 <- integrate(function(x) {
    ##   svt(ilogit(x), lower, upper, rho, nu) * ilogit(x) * (1-ilogit(x))
    ## }, lower = -100, upper = logit(.99999), subdivisions = 1000L)
    p1 <- integrate(function(x) bv_direct(x, upper[1], rho, nu),
                    lower = qfun(1e-5, nu),
                    upper = min(-qfun(1e-5, nu), upper[2]),
                    subdivisions = 1000L)
    p1 <- p1$value
    if(!has_df) {
      p2 <- mvtnorm::pmvnorm(lower, upper, corr = corr)
    } else {
      p2 <- mvtnorm::pmvt(lower, upper, df = nu, corr = corr)
    }
    p2 <- p2[1]
    ## p2 <- if(p2 < .5) log(p2) else log(1-p2)
    abs_err <- abs(p1 - p2)
    rel_err <- abs_err / abs(p1 + p2 + .1)
    expect_true(rel_err < 1e-3)
  }
}

curve(svt(x, lower, upper, rho, nu), from = 0, to = 1)

curve(svt(ilogit(x), lower, upper, rho, nu) * ilogit(x) * (1 - ilogit(x)),
      from = -100, to = 100)


#' # Alternative Calculation
#'
#' We are mainly interested in lower tails, which can be obtained via the one-dimensional integral
#' ```
#' int_{-inf}^b_1 Pr(X_2 < b_2 | X_1 = x) f(x) dx
#' ```

#' Direct CDF integrand for bivariate student-t.
#'
#' @param x Vector of integrand values.
#' @param b Upper quantile inside integrand
#' @param rho Correlation parameter.
#' @param nu Degrees of freedom parameter.  Setting to zero means normal CDF.
#'
#' @return Vector of same length as `x`.
bv_direct <- function(x, b, rho, nu) {
  loc <- rho * x
  scale <- 1 - rho^2
  if(nu > 0) {
    loc <- rho * x
    scale <- sqrt((nu + x^2)/(nu + 1) * scale)
    ans <- pt((b-loc)/scale, df = nu + 1) * dt(x, df = nu)
  } else {
    scale <- sqrt(scale)
    ans <- pnorm((b-loc)/scale) * dnorm(x)
  }
  ans
}

rho <- 0
nu <- runif(1, 0, 10)
lower <- c(-Inf, -Inf)
upper <- rnorm(2)

curve(svt(x, lower, upper, rho, nu), from = 0, to = 1)

curve(bv_direct(x, upper[1], rho, nu), from = -100, to = upper[2])

integrate(function(x) svt(x, c(-Inf, -Inf), upper, rho, nu),
          lower = 0, upper = 1, subdivisions = 1000L)

integrate(function(x) bv_direct(x, upper[1], rho, nu),
          lower = qt(1e-10, nu), upper = min(-qt(1e-10, nu), upper[2]),
          subdivisions = 1000L)
