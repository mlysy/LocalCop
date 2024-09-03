#--- unexported utility functions ----------------------------------------------

#' Construct the famil set.
#'
#' @details Calculates `nper` estimates of Kendall tau for non-overlaping sets of `u1` and `u2`.  If these are all positive, then can use any of the copula families.  Otherwise, can only use those which allow bot positive and negative dependence.
#'
#' **TODO:** Add rotated families for the case with only negative dependence.
#' @noRd
.get_family <- function(u1, u2, nper) {
  # empirical tau on non-overlapping periods
  etau <- .get_tau(u1, u2, nper)
  if(all(etau > 0)) {
    family  <- c(1:5, 13:14) # only +ve dependence
  } else if(all(etau < 0)) {
    family <- c(1,2,5, 23:24, 33:34) # only -ve dependence
  } else {
    family <- c(1,2,5) # +ve and -ve dependence
  }
  family
}

#' Kendall's tau on non-overlaping sets.
#'
#' @param u1, u2 Vectors of uniforms.
#' @param ntau Number of taus to calculate.
#' @return A vector of `ntau` taus.
#' @details Divides `u1` and `u2` into `ntau` non-overlaping sets and calculates the Kendall tau for each set.
#' @noRd
.get_tau <- function(u1, u2, ntau) {
  n <- length(u1)
  irng <- unique(round(seq(1, n, len = ntau+1)))
  irng <- cbind(irng[1:ntau], irng[-1])
  apply(irng, 1,
        function(rng) {
          ind <- rng[1]:rng[2]
          cor(u1[ind], u2[ind], method = "kendall")
        })
}

#' Get bandwidth set.
#'
#' @noRd
.get_band <- function(x, nband) {
  dx <- diff(sort(x),1)
  h.min <- max(dx)
  h.max <- max(x)-min(x)
  # get nband+2 values and remove smallest two
  log.seq <- seq(from=log(h.min), to=log(h.max), length.out = (nband+2))
  band <- round(exp(log.seq),5)
  band[-(1:2)]
}

#' Default optimization function.
#'
#' @noRd
.optim_default <- function(obj) {
  ## # coarse optimization: gradient-free
  ## opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
  ##              method = "Nelder-Mead",
  ##              control = list(maxit = 50, reltol = 1e-2))
  ## # fine optimization: quasi-newton (gradient-based)
  ## opt <- optim(par = opt$par, fn = obj$fn, gr = obj$gr,
  ##              method = "BFGS")
  opt <- stats::nlminb(start = obj$par,
                       objective = obj$fn,
                       gradient = obj$gr)
  # only need constant term since xc = 0 at x0 = x[ii]
  return(opt$par[1])
}

#' Estimate `eta` and/or `nu` if required.
#'
#' @param eta,nu Optional values of `eta` and/or `nu`.  If either of these is missing or `NA`, then uses [VineCopula::BiCopEst()] to estimate the parameters.
#' @noRd
.get_etaNu <- function(u1, u2, family, degree, eta, nu) {
  if(missing(eta)) eta <- NA
  if(missing(nu)) nu <- NA
  if(anyNA(eta) || (anyNA(nu) && family == 2)) {
    res <-  VineCopula::BiCopEst(u1 = u1, u2 = u2, family = family)
  }
  if(anyNA(eta)) {
    eta <- BiCopPar2Eta(family = family, par = res$par, par2 =res$par2)$eta
    if(degree == 1) eta <- c(eta, 0)
  }
  if(anyNA(nu)) {
    nu <- if(family == 2) res$par2 else 0
  }
  list(eta = eta, nu = nu)
}

#' Determine whether to run code in parallel.
#'
#' @noRd
.check_parallel <- function(cl) {
  pareval <- !anyNA(cl)
  if (!requireNamespace("parallel", quietly = TRUE)) {
    message("Package \"parallel\" needed for parallel evaluation.  Running serially instead.")
    pareval <- FALSE
  }
  pareval
}

#' Check that degree is valid.
#'
#' @noRd
.check_degree <- function(degree) {
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  ## degree <- match.arg(degree)
  ## return(as.numeric(degree == "linear"))
}

## #' Local likelihood estimation at a single covariate value.
## #'
## #' @inheritParams CondiCopLocFun
## #' @return List with elements \code{eta} and \code{loglik}.
## #' @export
## CondiCopLocFit1 <- function(u1, u2, family,
##                             X, x0, wgt, degree = c("linear", "constant"),
##                             eta, nu = 0) {
##   obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
##                         z = z, wgt = wgt, degree = degree, eta = eta, nu = nu)
##   opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
##                method = "BFGS")
##   return(list(eta = as.numeric(opt$par), loglik = -opt$value))
## }

#' Check whether copula family is known and/or supported.
#'
#' @noRd
.check_family <- function(family) {
  if(length(family) != 1) stop("family must be a single integer.  See `ConvertPar` for list of supported families.")
  if((family %in% .FamilySet) &&
     !(family %in% .FamilySet_supported)) {
    stop("Unsupported copula family.  See `?ConvertPar` for list of supported families.")
  } else if(!family %in% .FamilySet) {
    stop("Unknown copula family.  See `?ConvertPar` for list of supported families.")
  }
}
