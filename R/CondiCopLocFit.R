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
#' @param eta Optional initial value of the copula dependence parameter (scalar).  If missing will be estimated unconditionally by \code{VineCopula::BiCopEst}.
#' @param nu Optional initial value of other copula parameters (if they exist).  If missing will be estimated unconditionally by \code{VineCopula::BiCopEst}.
#' @param kernel Kernel function.  See \code{\link{KernFun}}.
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
  # initialize eta and nu
  degree <- match.arg(degree)
  etaNu <- .get_etaNu(u1 = u1, u2 = u2, family = family,
                      degree = degree, eta = eta, nu = nu)
  ieta <- etaNu$eta
  inu <- etaNu$nu
  fun <- function(xi) {
    wgt <- KernWeight(X = X, x = xi, band = band,
                      kernel = kernel, band.method = "constant")
    ## CondiCopLocFit1(u1 = u1, u2 = u2, family = family,
    ##                 z = X-xi, wgt = wgt, degree = degree,
    ##                 eta = ieta, nu = inu)$eta[1]
    obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
                          X = X, x = xi,
                          wgt = wgt, degree = degree, eta = ieta, nu = inu)
    # quasi-newton (gradient-based) optimization
    opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
                 method = "BFGS")
    return(as.numeric(opt$par[1])) # only need constant term since xc = 0 at x = xi
  }
  if(nx == 1) {
    eeta <- fun(x)
  } else {
    pareval <- .get_pareval(cl)
    if(!pareval) {
      # run serially
      eeta <- sapply(x, fun)
    } else {
      # run in parallel
      parallel::clusterExport(cl,
                              varlist = c("fun", "u1", "u2", "family", "X",
                                          "band", "kernel", "ieta", "inu"),
                              envir = environment())
      eeta <- parallel::parSapply(cl, X = x, FUN = fun)
    }
  }
  return(list(eta = eeta, nu = inu))
}

#--- helper functions ----------------------------------------------------------

# estimate eta and/or nu if required
.get_etaNu <- function(u1, u2, family, degree, eta, nu) {
  if(missing(eta) || (missing(nu) && family == 2)) {
    res <-  VineCopula::BiCopEst(u1 = u1, u2 = u2, family = family)
  }
  if(missing(eta)) {
    eta <- BiCopPar2Eta(family = family, par = res$par, par2 =res$par2)
    if(degree == "linear") eta <- c(eta, 0)
  }
  if(missing(nu)) {
    nu <- res$par2
  }
  list(eta = eta, nu = nu)
}

# determine whether to run code in parallel.
.get_pareval <- function(cl) {
  pareval <- !anyNA(cl)
  if (!requireNamespace("parallel", quietly = TRUE)) {
    message("Package \"parallel\" needed for parallel evaluation.  Running serially instead.")
    pareval <- FALSE
  }
  pareval
}
