#--- unexported utility functions ----------------------------------------------

# family set for given pair
# TODO: Add rotated families for the case with only negative dependence
.get_family <- function(u1, u2, nper) {
    # empirical tau on non-overlapping periods
    etau <- .get_tau(u1, u2, nper)
    if(all(etau > 0)) {
      family  <- 1:5   # case with only positive dependence
    } else {
      family <- c(1,2,5) # cases with both positive and negative dependence
    }
  }

# kendall's tau on non-overlaping sets
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

# get bandwidth set
.get_band <- function(X, nband) {
  dx <- diff(sort(X),1)
  h.min <- max(dx)
  h.max <- max(X)-min(X)
  # get nband+2 values and remove smallest two
  log.seq <- seq(from=log(h.min), to=log(h.max), length.out = (nband+2))
  band <- round(exp(log.seq),5)
  band[-(1:2)]
}

# default optimization function
.optim_default <- function(obj) {
  # coarse optimization: gradient-free
  opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
               method = "Nelder-Mead",
               control = list(maxit = 50, reltol = 1e-2))
  # fine optimization: quasi-newton (gradient-based)
  opt <- optim(par = opt$par, fn = obj$fn, gr = obj$gr,
               method = "BFGS")
  # only need constant term since xc = 0 at x = X[ii]
  return(opt$par[1])
}

# estimate eta and/or nu if required
# warning: passing NAs will result in fitting
.get_etaNu <- function(u1, u2, family, degree, eta, nu) {
  if(missing(eta)) eta <- NA
  if(missing(nu)) nu <- NA
  if(anyNA(eta) || (anyNA(nu) && family == 2)) {
    res <-  VineCopula::BiCopEst(u1 = u1, u2 = u2, family = family)
  }
  if(anyNA(eta)) {
    eta <- BiCopPar2Eta(family = family, par = res$par, par2 =res$par2)
    if(degree == "linear") eta <- c(eta, 0)
  }
  if(anyNA(nu)) {
    nu <- if(family == 2) res$par2 else 0
  }
  list(eta = eta, nu = nu)
}

# determine whether to run code in parallel.
.check_parallel <- function(cl) {
  pareval <- !anyNA(cl)
  if (!requireNamespace("parallel", quietly = TRUE)) {
    message("Package \"parallel\" needed for parallel evaluation.  Running serially instead.")
    pareval <- FALSE
  }
  pareval
}
