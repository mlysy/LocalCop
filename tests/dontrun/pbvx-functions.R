#' Calculate the integrand for rectangular BVX probabilities in algorithm QRSVT.
#'
#' Integrating this function gives the rectangular probabilities of a standard bivariate normal or student-t (BVX) distribution, i.e., which are marginally `N(0, 1)` or `t_(nu)`.
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
bvx_svt <- function(u, lower, upper, rho, nu) {
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
  .bvx_svt(u, lower = lims[1,], upper = lims[2,],
           rho = rho, nu = nu)
}

#' Raw calculation.
#'
#' Argument formatting is done in `svt()`.
#' @noRd
.bvx_svt <- function(u, lower, upper, rho, nu) {
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

#' Calculate the integrand for CDF of BVX using direct method.
#'
#' @param x Vector of integrand values.
#' @param b Upper quantile inside integrand.
#' @param rho Correlation parameter.
#' @param nu Degrees of freedom parameter.  Setting to zero means normal CDF.
#'
#' @return Vector of same length as `x`.
#'
#' @details Calculates the CDF via the integral
#' ```
#' int_{-inf}^b_1 Pr(X_2 < b_2 | X_1 = x) f(x) dx,
#' ```
#' where `(X_1, X_2)` are standard BVX variables and `f(x)` is the PDF of the standard normal or t distribution.
bvx_direct <- function(x, b, rho, nu) {
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

#' Integrand of BVN survival function using Drezner-Wesolowsky method.
#'
#' @param x Vector of integrand values.
#' @param b Upper quantile inside integrand.
#' @param rho Correlation parameter.
#'
#' @return Vector of the same length as `x`.
#' @details This is just the integrand.  An additional number is required to give the actual survival function.  The transformation for `|rho| ~ 1` not yet implemented.
