#' Local likelihood estimation.
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
#' @param optim_fun Optional optimization function to replace the default.  If provided, \code{optim_fun} should take a single argument corresponding the output of \code{\link{CondiCopLocFun}}, and return a scalar value corresponding to the estimate of \code{eta} at a given covariate value in \code{x}.
#' @param cl Optional parallel cluster created with \code{parallel::makeCluster}.  If \code{NA} runs serially.
#' @return Vector of estimated dependence parameter \code{eta} at each value of \code{x}.
#' @export
CondiCopLocFit <- function(u1, u2, family, X, x, nx = 100,
                           degree = c("linear", "constant"),
                           eta, nu, kernel = KernEpa, band,
                           optim_fun, cl = NA) {
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
  # optimization function
  if(missing(optim_fun)) {
    optim_fun <- .optim_default
  }
  fun <- function(xi) {
    wgt <- KernWeight(X = X, x = xi, band = band,
                      kernel = kernel, band.method = "constant")
    obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
                          X = X, x = xi,
                          wgt = wgt, degree = degree, eta = ieta, nu = inu)
    return(optim_fun(obj))
  }
  if(nx == 1) {
    eeta <- fun(x)
  } else {
    if(!.check_parallel(cl)) {
      # run serially
      eeta <- sapply(x, fun)
    } else {
      # run in parallel
      parallel::clusterExport(cl,
                              varlist = c("fun", "u1", "u2", "family", "X",
                                          "band", "kernel", "optim_fun",
                                          "ieta", "inu"),
                              envir = environment())
      eeta <- parallel::parSapply(cl, X = x, FUN = fun)
    }
  }
  return(list(eta = as.numeric(eeta), nu = as.numeric(inu)))
}

#--- helper functions ----------------------------------------------------------

