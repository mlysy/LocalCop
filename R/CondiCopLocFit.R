#' Local likelihood estimation.
#'
#' Estimate the bivariate copula dependence parameter `eta` at multiple covariate values.
#'
#' @template param-u1
#' @template param-u2
#' @template param-family
#' @template param-X
#' @template param-xseq
#' @param nx If `x` is missing, defaults to `nx` equally spaced values in `range(X)`.
#' @template param-degree
#' @param eta Optional initial value of the copula dependence parameter (scalar).  If missing will be estimated unconditionally by [VineCopula::BiCopEst()].
#' @param nu Optional initial value of second copula parameter, if it exists.  If missing and required, will be estimated unconditionally by [VineCopula::BiCopEst()].  If provided and required, will not be estimated.
#' @template param-kernel
#' @template param-band
#' @param optim_fun Optional specification of local likelihood optimization algorithm.  See **Details**.
#' @param cl Optional parallel cluster created with [parallel::makeCluster()], in which case optimization for each element of `x` will be done in parallel on separate cores.  If `cl == NA`, computations are run serially.
#' @return List with the following elements:
#' \describe{
#'   \item{`x`}{The vector of covariate values at which the local likelihood is fit.}
#'   \item{`eta`}{The vector of estimated dependence parameters of the same length as `x`.}
#'   \item{`nu`}{The scalar value of the estimated (or provided) second copula parameter.}
#' }
#' @details By default, optimization is performed by taking a few steps with the gradient-free simplex algorithm `stats::optim(method = "Nelder-Mead")` for stability, then continuing with the quasi-Newton algorithm `stats::optim(method = "BFGS")`, which uses gradient information provided by automatic differentiation (AD) as implemented by \pkg{TMB}.
#'
#' If the default method is to be overridden, `optim_fun` should be provided as a function taking a single argument corresponding to the output of [CondiCopLocFun()], and return a scalar value corresponding to the estimate of `eta` at a given covariate value in `x`.  Note that \pkg{TMB} calculates the *negative* local (log)likelihood, such that the objective function is to be minimized.  See **Examples**.
#' @example examples/CondiCopLocFit.R
#' @export
CondiCopLocFit <- function(u1, u2, family, X, x, nx = 100,
                           degree = 1,
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
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  ## degree <- match.arg(degree)
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
                      kernel = kernel, band_type = "constant")
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
  return(list(x = x, eta = as.numeric(eeta), nu = as.numeric(inu)))
}

