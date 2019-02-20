#' Local likelihood estimation at a single covariate value.
#'
#' Optionally performs bandwidth and family selection.
#'
#' @param u1 Vector of first uniform response.
#' @param u2 Vector of second uniform response.
#' @template param-family
#' @param X Vector of observed covariate values.
#' @param x Vector of covariate values within \code{range(X)} at which to fit the local likelihood.  Does not have to be a subset of \code{X}.
#' @param nx If \code{x} is missing, defaults to \code{nx} equally spaced values in \code{range(X)}.
#' @param degree Character vector specifying degree of local polynomial.  Either "constant" or "linear".
#' @param eta Optional initial value of the copula dependence parameter (scalar).
#' @param nu Optional initial value of other copula parameters (scalar or vector).
#' @param kernel Kernel function.  See \code{\link{kernel}}.
#' @param band Bandwidth parameter (positive scalar).
#' @param cl Optional parallel cluster created with \code{parallel::makeCluster}.  If \code{NA} runs serially.
#' @return Vector of estimated dependence parameter \code{eta} at each value of \code{x}.
#' @note FIXME: \code{parallel} requirement.
#' @export
CondiCopLocFit <- function(u1, u2, family, X, x, nx = 100,
                           degree = c("linear", "constant"),
                           eta, nu,
                           kernel = KernEpa, band, cl = NA) {
  # default x
  if(missing(x)) {
    x <- seq(min(X), max(X), len = nx)
  } else {
    x <- sort(x)
  }
  nx <- length(x)
  # default eta and nu
  degree <- match.arg(degree)
  if(missing(eta) || missing(nu)) {
    res <-  VineCopula::BiCopEst(u1,u2,family)
    if(missing(eta)) {
      eta <- BiCopPar2Eta(family=family, par = res$par, par2 =res$par2)
      if(degree == "linear") eta <- c(eta, 0)
    }
    if(missing(nu)) {
      nu <- res$par2
    }
  }
  # detect parallel
  fun <- function(xi) {
    wgt <- KernWeight(X = X, x = xi, band = band,
                      kernel = kernel, band.method = "constant")
    CondiCopLocFit1(u1 = u1, u2 = u2, family = family,
                    z = X-xi, wgt = wgt, degree = degree,
                    eta = eta, nu = nu)$eta[1]
  }
  if(anyNA(cl)) {
    eta.hat <- sapply(x, fun)
  } else {
    clusterExport(cl,
                  varlist = c("fun", "u1", "u2", "family", "X",
                              "band", "kernel", "eta", "nu"),
                  envir = environment())
    eta.hat <- parSapply(cl, X = x, FUN = fun)
  }
  return(list(eta = eta.hat, nu = nu))
}
