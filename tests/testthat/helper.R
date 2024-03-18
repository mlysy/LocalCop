set.seed(1772)

#' Generate data and parameter values for local likelihood.
#'
#' @param family Copula family: integer 1-5.
#' @return List with elements `udata`, `epar`, `epar2`, `wgt`, `x`, `x0`, `eta`.
data_sim <- function(family) {
  # generate data
  nobs <- sample(10:50, 1)  # number of observations
  x <- sort(runif(nobs))    # values in  (0,1] interval
  # etafun <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
  etafun <- function(t) 2*t   # linear calibration function
  tpar <- BiCopEta2Par(family = family, etafun(x))$par # true copula parameter
  tpar2 <- 10 + runif(1) # only relevant for two-parameter copulas
  udata <- VineCopula::BiCopSim(N=nobs, family=family,
                                par = tpar, par2 = tpar2)
  # local likelihood
  x0 <- runif(1, min(x), max(x)) # evaluation point
  # weight specification
  kern <- sample(c(KernEpa, KernGaus, KernBeta,
                   KernBiQuad, KernTriAng), 1)[[1]]
  band <- runif(1, .025, .5)
  wgt <- KernWeight(x = x, x0 = x0, band = band, kernel = kern)
  ## wgt <- rep(1, nobs)
  # local likelihood calculation
  for(ii in 1:100) {
    # generate valid eta/epar pair:
    # for family == 3 need epar < 27.9
    # for family == 4 need epar < 16.9
    eta <- rnorm(2)/2  # evaluation parameter
    epar <- BiCopEta2Par(family = family, eta = eta[1] + eta[2] * (x-x0))$par
    if((!family %in% c("3", "4")) ||
       ((family == "3") && (max(epar) < 27.9)) ||
       ((family == "4") && (max(epar) < 16.9))) {
      break
    }
  }
  ## eta <- rnorm(2)/2  # evaluation parameter
  ## epar <- BiCopEta2Par(family = family, eta = eta[1] + eta[2] * (x-x0))$par
  ## if(family == "3") {
  ##   # restricted range in VineCopula package
  ##   epar <- pmin(epar, 27.9)
  ## }
  ## if(family == "4") {
  ##   # restricted range in VineCopula package
  ##   epar <- pmin(epar, 16.9)
  ## }
  epar2 <- 10 + runif(nobs) # only relevant for two parameter copulas
  list(udata = udata,
       epar = epar, epar2 = epar2, wgt = wgt,
       x = x, x0 = x0, eta = eta)
}
