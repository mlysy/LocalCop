#' Bandwidth and/or family selection.
#'
#' Selects among a set of bandwidths and/or copula families the one which maximizes the cross-validated local likelihood.
#'
#' @inheritParams CondiCopLikCV
#' @param xind Specification of \code{xind} for each bandwidth.  Can be a scalar integer, a vector of \code{nband} integers, or a list of \code{nband} vectors of integers.
#' @param nu Optional vector of fixed \code{nu} parameter for each family.  If missing or \code{NA} get estimated from the data.
#' @param nband Number of pilot bandwidth values used in selection.
#' @return A list with elements \code{band} and \code{family} with their optimal values.
#' @export
CondiCopSelect <- function(u1, u2, family, X, xind = 100,
                           degree = c("linear", "constant"),
                           nu,
                           kernel = KernEpa, band, nband = 6,
                           optim_fun, cl = NA) {
  ## degree <- match.arg(degree)
  ## n <- length(u1)
  # family set
  if(missing(family)) family <- .get_family(u1, u2, nper = 10)
  nfam <- length(family)
  # initial parameters
  if(missing(nu)) nu <- rep(NA, nfam)
  nu <- sapply(1:nfam, function(ii) {
    .get_etaNu(u1 = u1, u2 = u2, family = family[ii],
               degree = degree, eta = c(1,0), nu = nu[ii])$nu
  })
  # bandwidth set
  if(missing(band)) band <- .get_band(X, nband)
  nband <- length(band)
  # selection process
  if(nband == 1 & length(family)==1) {
    # no need to do selection if there is only one choice.
    res <- list(band=band, family=family)
  } else {
    # calculate cvLIK for family & bandwidth combinations
    gridVal <- expand.grid(band = band, family = family)
    gridVal <- cbind(gridVal, nu = rep(nu, each = nband))
    # xind value adjustment for CV calculation
    if(is.numeric(xind)) {
      if(length(xind) == 1) xind <- rep(xind, nband)
      xind <- as.list(xind)
    }
    xind <- rep(xind, nfam)
    if(length(xind) != nrow(gridVal)) {
      stop("Incorrect specification of xind.")
    }
    # optimization function
    if(missing(optim_fun)) {
      optim_fun <- .optim_default
    }
    fun <- function(ii) {
      CondiCopLikCV(u1=u1, u2=u2, family = gridVal$family[ii],
                    X=X, xind = xind[[ii]], degree = degree,
                    eta=c(1,0), nu=nu[ii], kernel=kernel,
                    band = gridVal$band[ii],
                    optim_fun = optim_fun, cl = NA)
    }
    if(!.check_parallel(cl)) {
      # run serially
      cvLIK <- sapply(1:nrow(gridVal), fun)
    } else {
      # run in parallel
      parallel::clusterExport(cl,
                              varlist = c("fun", "u1", "u2", "family", "X",
                                          "band", "kernel", "optim_fun",
                                          "gridVal", "xind"),
                              envir = environment())
      cvLIK <- parallel::parSapply(cl, X = 1:nrow(gridVal), FUN = fun)
    }
    isel <- which.max(cvLIK)
    res <- list(family = gridVal$family[isel],
                band = gridVal$band[isel])
  }
  return(res)
}

#--- helper functions -----------------------------------------------------

