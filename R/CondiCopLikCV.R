#' Cross-validated likelihood.
#'
#' Leave-one-out local likelihood copula parameter estimates are interpolated, then used to calculate the conditional copula likelihood function.
#'
#' @template param-u1
#' @template param-u2
#' @template param-family
#' @template param-X
#' @param xind Vector of indices in \code{X} at which to calculate leave-one-out parameter estimates.  Can also be supplied as a single integer, in which case \code{xind} equally spaced observations are taken from \code{X}.
#' @template param-degree
#' @param eta,nu,kernel,band,optim_fun,cl See \code{\link{CondiCopLocFit}}.
#' @param cveta_out If \code{TRUE}, the CV estimate of eta at each point in \code{X} in addition to the CV log-likelihood.
#' @return Scalar value of the cross-validated log-likelihood, or a list with elements \code{eta}, \code{nu} and \code{loglik} if \code{cveta_out = TRUE}.
#' @example examples/CondiCopLikCV.R
#' @export
CondiCopLikCV <- function(u1, u2, family, X, xind = 100,
                          degree = 1,
                          eta, nu, kernel = KernEpa, band,
                          optim_fun, cveta_out = FALSE, cl = NA) {
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
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  ## degree <- match.arg(degree)
  etaNu <- .get_etaNu(u1 = u1, u2 = u2, family = family,
                      degree = degree, eta = eta, nu = nu)
  ieta <- etaNu$eta
  inu <- etaNu$nu
  # cross validation: estimation step
  # optimization function
  if(missing(optim_fun)) {
    optim_fun <- .optim_default
  }
  fun <- function(ii) {
    wgt <- KernWeight(X = X[-ii], x = X[ii], band = band,
                      kernel = kernel, band_type = "constant")
    obj <- CondiCopLocFun(u1 = u1[-ii], u2 = u2[-ii], family = family,
                          X = X[-ii], x = X[ii],
                          wgt = wgt, degree = degree, eta = ieta, nu = inu)
    return(optim_fun(obj))
  }
  if(!.check_parallel(cl)) {
    # run serially
    cveta <- sapply(xind, fun)
  } else {
    # run in parallel
    parallel::clusterExport(cl,
                            varlist = c("fun", "u1", "u2", "family", "X",
                                        "band", "kernel", "optim_fun",
                                        "ieta", "inu"),
                            envir = environment())
    cveta <- parallel::parSapply(cl, X = xind, FUN = fun)
  }
  # validation step
  cveta <- approx(X[xind],
                  y = cveta, xout = X)$y # interpolate cveta to all observations
  obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
                        X = cveta, x = 0, eta = c(0,1), nu = inu,
                        wgt = rep(1, length(u1)), degree = 1)
  cvll <- -obj$fn(c(0,1))
  # correct for likelihood constants
  if(family == 2) {
    # Student-t
    cst <- lgamma(.5*(inu+2)) + lgamma(.5*inu) - 2*lgamma(.5*(inu+1))
    cvll <- cvll + length(X) * cst
  }
  if(family == 4) {
    cvll <- cvll - sum(log(u1) + log(u2))
  }
  if(!cveta_out) {
    return(cvll)
  } else {
    return(list(eta = cveta, nu = inu, loglik = cvll))
  }
}

#--- scratch -------------------------------------------------------------------

## plot_fun <- function(eta0, eta1, npts = 100) {
##   eta0 <- seq(eta0[1], eta0[2], len = npts)
##   eta1 <- seq(eta1[1], eta1[2], len = npts)
##   Eta <- as.matrix(expand.grid(eta0, eta1))
##   ld <- -apply(Eta, 1, obj$fn)
##   ld <- matrix(ld-range(ld,na.rm=TRUE,finite=TRUE)[2], npts, npts)
##   ## ld[is.na(ld)] <- min(abs(ld[!is.na(ld)])) + 1e5
##   contour(x = eta0, y = eta1, z = exp(ld))
## }
