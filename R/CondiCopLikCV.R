#' Cross-validated likelihood.
#'
#' Leave-one-out local likelihood copula parameter estimates are interpolated, then used to calculate the conditional copula likelihood function.
#'
#' @template param-u1
#' @template param-u2
#' @template param-family
#' @template param-x
#' @param xind Vector of indices in `sort(x)` at which to calculate leave-one-out parameter estimates.  Can also be supplied as a single integer, in which case `xind` equally spaced observations are taken from `x`.
#' @template param-degree
#' @param eta,nu,kernel,band,optim_fun,cl See [CondiCopLocFit()].
#' @template param-cv_all
#' @param cveta_out If `TRUE`, return the CV estimate of eta at each point in `x` in addition to the CV log-likelihood.
#' @return If `cveta_out = FALSE`, scalar value of the cross-validated log-likelihood.  Otherwise, a list with elements:
#' \describe{
#'   \item{`x`}{The sorted values of `x`.}
#'   \item{`eta`}{The leave-one-out estimates interpolated from the values in `xind` to all of those in `x`.}
#'   \item{`nu`}{The scalar value of the estimated (or provided) second copula parameter.}
#'   \item{`loglik`}{The cross-validated log-likelihood.}
#' }
#' @seealso This function is typically used in conjunction with [CondiCopSelect()]; see example there.
#' @export
CondiCopLikCV <- function(u1, u2, family, x, xind = 100,
                          degree = 1,
                          eta, nu, kernel = KernEpa, band,
                          optim_fun, cveta_out = FALSE,
                          cv_all = FALSE, cl = NA) {
  # sort observations
  ix <- order(x)
  x <- x[ix]
  u1 <- u1[ix]
  u2 <- u2[ix]
  # index of validation observations
  if(length(xind) == 1) {
    xind <- unique(round(seq(1, length(x), len = xind)))
  }
  # initialize eta and nu
  .check_family(family)
  .check_degree(degree)
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
    wgt <- KernWeight(x = x[-ii], x0 = x[ii], band = band,
                      kernel = kernel, band_type = "constant")
    obj <- CondiCopLocFun(u1 = u1[-ii], u2 = u2[-ii], family = family,
                          x = x[-ii], x0 = x[ii],
                          wgt = wgt, degree = degree, eta = ieta, nu = inu)
    return(optim_fun(obj))
  }
  if(!.check_parallel(cl)) {
    # run serially
    cveta <- sapply(xind, fun)
  } else {
    # run in parallel
    parallel::clusterExport(cl,
                            varlist = c("fun", "u1", "u2", "family", "x",
                                        "band", "kernel", "optim_fun",
                                        "ieta", "inu"),
                            envir = environment())
    cveta <- parallel::parSapply(cl, X = xind, FUN = fun)
  }
  # validation step
  # interpolate cveta to all observations
  cveta <- approx(x[xind], y = cveta, xout = x)$y
  if(cv_all) xind <- 1:length(u1)
  nx <- length(xind)
  obj <- CondiCopLocFun(u1 = u1[xind], u2 = u2[xind], family = family,
                        x = cveta[xind], x0 = 0, eta = c(0,1), nu = inu,
                        wgt = rep(1, nx), degree = 1)
  cvll <- -obj$fn(c(0,1))
  ## # correct for likelihood constants
  ## if(family == 2) {
  ##   # Student-t
  ##   cst <- lgamma(.5*(inu+2)) + lgamma(.5*inu) - 2*lgamma(.5*(inu+1))
  ##   cvll <- cvll + nx * cst
  ## }
  ## if(family == 4) {
  ##   cvll <- cvll - sum(log(u1[xind]) + log(u2[xind]))
  ## }
  if(!cveta_out) {
    return(cvll)
  } else {
    return(list(x = x, eta = cveta, nu = inu, loglik = cvll))
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
