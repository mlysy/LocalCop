#--- test local likelihood implementation in TMB -------------------------------

context("LocLikFun")

test_that("LocLikFun is same in VineCopula and TMB", {
  nreps <- 100
  for(ii in 1:nreps) {
    # generate data
    n <- sample(10:50, 1) # number of time points
    X <- sort(runif(n))  # time points in  (0,1] interval
    family <- sample(1:5, 1)
    eta.fnc <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
    tpar <- BiCopEta2Par(family = family, eta.fnc(X))$par
    tpar2 <- 10 + runif(1)
    udata <- VineCopula::BiCopSim(N=n, family=family,
                                  par = tpar, par2=tpar2)
    # local likelihood
    x0 <- runif(1, min(X), max(X)) # evaluation point
    # weight specification
    kern <- sample(c(KernEpa, KernGaus, KernBeta,
                     KernBiQuad, KernTriAng), 1)[[1]]
    band <- runif(1, .025, .5)
    wgt <- KernWeight(X = X, x = x0, band = band, kernel = kern)
    # local likelihood calculation
    eeta <- rnorm(2)/2
    ## eta <- eeta[1] + eeta[2] * (X-x0)
    ## if(family == 4) {
    ##   # gumbel parameter must be between [1, 17]
    ##   eta <- pmin(pmax(eta, 1), 17)
    ## }
    epar <- BiCopEta2Par(family = family, eta = eeta[1] + eeta[2] * (X-x0))$par
    epar2 <- 10 + runif(1)
    # in R
    ll_r <- VineCopula::BiCopPDF(u1 = udata[,1], u2 = udata[,2],
                                 family = family, par = epar, par2 = epar2)
    ll_r <- sum(wgt * log(ll_r))
    # in TMB
    ll_tmb <- CondiCopLocFun(u1 = udata[,1], u2 = udata[,2],
                             family = family, X = X, x = x0,
                             wgt = wgt, eta = eeta, nu = epar2)
    ll_tmb <- -ll_tmb$fn(eeta)
    if(family == 2) {
      # didn't include factor of nu in TMB as we don't optimize wrt it
      cst <- lgamma(.5*(epar2+2)) + lgamma(.5*epar2) - 2*lgamma(.5*(epar2+1))
      ll_tmb <-  ll_tmb + sum(cst * wgt)
    }
    if(family == 4) {
      # extra factor of log(u1) and log(u2)
      cst <- -(log(udata[,1]) + log(udata[,2]))
      ll_tmb <- ll_tmb + sum(cst * wgt)
    }
    expect_equal(ll_r, ll_tmb)
  }
})
