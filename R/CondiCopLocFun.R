#' Create a \pkg{TMB} local likelihood function.
#'
#' Wraps a call to [TMB::MakeADFun()].
#'
#' @template param-u1
#' @template param-u2
#' @template param-family
#' @template param-x
#' @param x0 Scalar covariate value at which to evaluate the local likelihood.  Does not have to be a subset of `x`.
#' @param wgt Vector of positive kernel weights.
#' @template param-degree
#' @param eta Value of the copula dependence parameter.  Scalar or vector of length two, depending on whether `degree` is 0 or 1.
#' @param nu Value of the other copula parameter.  Scalar or vector of same length as `u1`.  Ignored if `family != 2`.
#' @return A list as returned by a call to [TMB::MakeADFun()].  In particular, this contains elements `fun` and `gr` for the *negative* local likelihood and its gradient with respect to `eta`.
#' @example examples/CondiCopLocFun.R
#' @export
CondiCopLocFun <- function(u1, u2, family,
                           x, x0, wgt, degree = 1,
                           eta, nu) {
  if(!family %in% 1:5) {
    stop("Unsupported copula family (must be integer between 1-5).")
  }
  wpos <- wgt > 0 # index of positive weights
  # create TMB function
  # format nu
  if(family != 2) nu <- 0 # second copula parameter
  if(length(nu) == 1) nu <- rep(nu, length(wgt))
  if(length(nu) != length(wgt)) {
    stop("nu must be of length 1 or have same length as wgt.")
  }
  # data input
  data <- list(model = "LocalLikelihood",
               y1 = u1[wpos], y2 = u2[wpos],
               wgt = wgt[wpos], xc = x[wpos]-x0,
               family = family, nu = nu[wpos])
  parameters <- list(beta = eta)
  # convert degree to TMB::map
  ## degree <- .format_degree(degree)
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  map <- list(beta = factor(c(1, 2)))
  if(degree == 0) {
    parameters$beta[2] <- 0
    map$beta[2] <- NA
  }
  TMB::MakeADFun(
    data = data,
    parameters = parameters,
    map = map,
    DLL = "LocalCop_TMBExports",
    silent = TRUE
  )
}


