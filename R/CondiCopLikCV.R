#' Cross-validated likelihood.
#'
#' @inheritParams CondiCopLocFit
#' @param xind Vector of indices in \code{X} at which to calculate leave-one-out parameter estimates.  Can also be supplied as a single integer, in which case \code{xind} equally spaced observations are taken from \code{X}.
#' @return Scalar value of the cross-validated log-likelihood.
#' @details Explain exactly what this does...
#' @export
CondiCopLikCV <- function(u1, u2, family, X, xind = 100,
                       degree = c("linear", "constant"),
                       eta, nu,
                       kernel = KernEpa, band, cl = NA) {
  # sort observations
  ix <- order(X)
  X <- X[ix]
  u1 <- u1[ix]
  u2 <- u2[ix]
  # index of validation observations
  if(length(xind) == 1) {
    xind <- unique(round(seq(1, length(X), len = xind)))
  }
  nx <- length(xind)
  # initialize eta and nu
  degree <- match.arg(degree)
  etaNu <- .get_etaNu(u1 = u1, u2 = u2, family = family,
                      degree = degree, eta = eta, nu = nu)
  ieta <- etaNu$eta
  inu <- etaNu$nu
  # cross validation: estimation step
  fun <- function(ii) {
    wgt <- KernWeight(X = X[-ii], x = X[ii], band = band,
                      kernel = kernel, band.method = "constant")
    ## CondiCopLocFit1(u1 = u1[-ii], u2 = u2[-ii], family = family,
    ##                 z = X[-ii]-X[ii], wgt = wgt, degree = degree,
    ##                 eta = eta, nu = nu)$eta[1]
    obj <- CondiCopLocFun(u1 = u1[-ii], u2 = u2[-ii], family = family,
                          X = X[-ii], x = X[ii],
                          wgt = wgt, degree = degree, eta = ieta, nu = inu)
    # quasi-newton (gradient-based) optimization
    opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
                 method = "BFGS")
    return(as.numeric(opt$par[1])) # only need constant term since xc = 0 at x = X[ii]
  }
  pareval <- .get_pareval(cl) # check if we can run in parallel
  if(!pareval) {
    # run serially
    cveta <- sapply(xind, fun)
  } else {
    # run in parallel
    parallel::clusterExport(cl,
                            varlist = c("fun", "u1", "u2", "family", "X",
                                        "band", "kernel", "ieta", "inu"),
                            envir = environment())
    cveta <- parallel::parSapply(cl, X = xind, FUN = fun)
  }
  # validation step
  cveta <- approx(X[xind],
                  y = cveta, xout = X)$y # interpolate cveta to all observations
  obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
                        X = cveta, x = 0, eta = c(0,1), nu = inu,
                        wgt = rep(1, length(u1)), degree = "linear")
  return(-obj$fn(c(0,1)))
}
