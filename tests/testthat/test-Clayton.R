#--- Clayton copula tests ------------------------------------------------------

##############
# PDF TEST
##############
test_that("ClaytonNLL is same in VineCopula and TMB", {
  nreps <- 100
  for(ii in 1:nreps) {
    # generate data
    nobs <- sample(10:50, 1) # number of observations
    X <- sort(runif(nobs))  # values in  (0,1] interval
    ## family <- sample(1:5, 1)
    family <- 3 # Clayton copula
    eta_fun <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
    tpar <- BiCopEta2Par(family = family, eta_fun(X))$par # true copula parameter
    tpar2 <- 10 + runif(1) # not relevant for Clayton
    udata <- VineCopula::BiCopSim(N=nobs, family=family,
                                  par = tpar, par2=tpar2)
    # local likelihood
    x0 <- runif(1, min(X), max(X)) # evaluation point
    # weight specification
    kern <- sample(c(KernEpa, KernGaus, KernBeta,
                     KernBiQuad, KernTriAng), 1)[[1]]
    band <- runif(1, .025, .5)
    wgt <- KernWeight(X = X, x = x0, band = band, kernel = kern)
    ## wgt <- rep(1, nobs)
    # local likelihood calculation
    eeta <- rnorm(2) # evaluation parameter
    epar <- BiCopEta2Par(family = family, eta = eeta[1] + eeta[2] * (X-x0))$par
    epar2 <- 10 + runif(1) # not relevant for Clayton
    # in R
    ll_r <- VineCopula::BiCopPDF(u1 = udata[,1], u2 = udata[,2],
                                 family = family, par = epar, par2 = epar2)
    ll_r <- -sum(wgt * log(ll_r))
    # in TMB
    cop_adf <- TMB::MakeADFun(
           data = list(
             model = "ClaytonNLL",
             u1 = udata[,1],
             u2 = udata[,2],
             weights = wgt
            ),
           parameters = list(theta = epar),
           silent = TRUE, DLL = "LocalCop_TMBExports")
    ll_tmb <- cop_adf$fn(epar)
    expect_equal(ll_r, ll_tmb)
  }
})


##############
# CDF TEST
##############
test_that("ClaytonCDF is same in VineCopula and TMB", {
  nreps <- 100
  for(ii in 1:nreps) {
    # generate data
    nobs <- sample(10:50, 1) # number of observations
    X <- sort(runif(nobs))  # values in  (0,1] interval
    ## family <- sample(1:5, 1)
    family <- 3 # Clayton copula
    eta_fun <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
    tpar <- BiCopEta2Par(family = family, eta_fun(X))$par # true copula parameter
    tpar2 <- 10 + runif(1) # not relevant for Clayton
    udata <- VineCopula::BiCopSim(N=nobs, family=family,
                                  par = tpar, par2=tpar2)
    # local likelihood
    x0 <- runif(1, min(X), max(X)) # evaluation point
    # weight specification
    kern <- sample(c(KernEpa, KernGaus, KernBeta,
                     KernBiQuad, KernTriAng), 1)[[1]]
    band <- runif(1, .025, .5)
    wgt <- KernWeight(X = X, x = x0, band = band, kernel = kern)
    ## wgt <- rep(1, nobs)
    # local likelihood calculation
    eeta <- rnorm(2) # evaluation parameter
    epar <- BiCopEta2Par(family = family, eta = eeta[1] + eeta[2] * (X-x0))$par
    epar2 <- 10 + runif(1) # not relevant for Clayton
    # in R
    ll_r <- VineCopula::BiCopCDF(u1 = udata[,1], u2 = udata[,2],
                                 family = family, par = epar, par2 = epar2)
    ll_r <- -sum(wgt * log(ll_r))
    # in TMB
    cop_adf <- TMB::MakeADFun(
      data = list(
        model = "ClaytonCDF",
        u1 = udata[,1],
        u2 = udata[,2],
        weights = wgt
      ),
      parameters = list(theta = epar),
      silent = TRUE, DLL = "LocalCop_TMBExports")
    ll_tmb <- cop_adf$fn(epar)
    expect_equal(ll_r, ll_tmb)
  }
})



