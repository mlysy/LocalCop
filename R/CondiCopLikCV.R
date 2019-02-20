#' Cross-validated likelihood.
#'
#' @inheritParams CondiCopLocFit
#' @param xind Vector of indices in \code{X} at which to calculate leave-one-out parameter estimates.  Can also be supplied as a single integer, in which case \code{xind} equally spaced observations are taken from \code{X}.
#' @return Scalar value of the cross-validated log-likelihood.
#' @details Explain exactly what this does...
#' @export
CondiCopLikCV <- function(u1, u2, family, X, xind,
                       degree = c("linear", "constant"),
                       eta, nu,
                       kernel = KernEpa, band, cl = NA) {
  # sort observations
  ix <- order(X)
  X <- X[ix]
  u1 <- u1[ix]
  u2 <- u2[ix]
  # index of validation observations
  if(missing(xind)) {
    xind <- unique(round(seq(1, length(X), len = nx)))
  }
  nx <- length(xind)
  # default eta and nu
  degree <- match.arg(degree)
  en <- .get_etaNu(u1, u2, family, degree, eta, nu)
  eta <- en$eta
  nu <- en$nu
  # cross validation: estimation step
  fun <- function(ii) {
    wgt <- KernWeight(X = X[-ii], x = X[ii], band = band,
                      kernel = kernel, band.method = "constant")
    CondiCopLocFit1(u1 = u1[-ii], u2 = u2[-ii], family = family,
                    z = X[-ii]-X[ii], wgt = wgt, degree = degree,
                    eta = eta, nu = nu)$eta[1]
  }
  if(anyNA(cl)) {
    cveta <- sapply(xind, fun)
  } else {
    clusterExport(cl,
                  varlist = c("fun", "u1", "u2", "family", "X",
                              "band", "kernel", "eta", "nu"),
                  envir = environment())
    cveta <- parSapply(cl, X = xind, FUN = fun)
  }
  # validation step
  # interpolate cveta to all observations
  cveta <- approx(X[xind], cveta, xout=X)$y
  obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
                        z = cveta, eta = c(0,1), nu = nu,
                        wgt = rep(1, length(u1)), degree = "linear")
  return(-obj$fn(c(0,1)))
}
