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
    ## p1 <- integrate(function(x) bvx_svt(x, lower, upper, rho, nu),
    ##                 lower = 0, upper = 1, subdivisions = 1000L)
    ## p1 <- integrate(function(x) {
    ##   bvx_svt(ilogit(x), lower, upper, rho, nu) * ilogit(x) * (1-ilogit(x))
    ## }, lower = -100, upper = logit(.99999), subdivisions = 1000L)
    p1 <- integrate(function(x) bvx_direct(x, upper[1], rho, nu),
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

curve(bvx_svt(x, lower, upper, rho, nu), from = 0, to = 1)

curve(bvx_svt(ilogit(x), lower, upper, rho, nu) * ilogit(x) * (1 - ilogit(x)),
      from = -100, to = 100)


#' # Alternative Calculation
#'
#' We are mainly interested in lower tails, which can be obtained via the one-dimensional integral
#' ```
#' int_{-inf}^b_1 Pr(X_2 < b_2 | X_1 = x) f(x) dx
#' ```

rho <- 0
nu <- runif(1, 0, 10)
lower <- c(-Inf, -Inf)
upper <- rnorm(2)

curve(bvx_svt(x, lower, upper, rho, nu), from = 0, to = 1)

curve(bvx_direct(x, upper[1], rho, nu), from = -100, to = upper[2])

integrate(function(x) bvx_svt(x, c(-Inf, -Inf), upper, rho, nu),
          lower = 0, upper = 1, subdivisions = 1000L)

integrate(function(x) bvx_direct(x, upper[1], rho, nu),
          lower = qt(1e-10, nu), upper = min(-qt(1e-10, nu), upper[2]),
          subdivisions = 1000L)
