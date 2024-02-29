#--- Copula tests ------------------------------------------------------

data_sim <- function(family) {
  # generate data
  nobs <- sample(10:50, 1)  # number of observations
  X <- sort(runif(nobs))    # values in  (0,1] interval
  # etafun <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
  etafun <- function(t) 2*t   # linear calibration function
  tpar <- BiCopEta2Par(family = family, etafun(X))$par # true copula parameter
  tpar2 <- 10 + runif(1) # only relevant for two-parameter copulas
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
  eeta <- rnorm(2)  # evaluation parameter
  epar <- BiCopEta2Par(family = family, eta = eeta[1] + eeta[2] * (X-x0))$par
  if(family == "3"){epar <- pmin(epar, 27.9)} # restricted range in VineCopula package
  if(family == "4"){epar <- pmin(epar, 16.9)} # restricted range in VineCopula package
  epar2 <- 10 + runif(1) # only relevant for two parameter copulas
  list(udata = udata, epar = epar, epar2 = epar2, wgt = wgt)
}



##############
# PDF TEST
##############

test_that("Copula density is same in VineCopula and TMB", {
  nreps <- 20
  test_descr <- expand.grid(family = c(3, 4, 5), # add copula families 
                            stringsAsFactors = FALSE)
  n_test <- nrow(test_descr)
  for(ii in 1:n_test) {
    for(jj in 1:nreps) {
      # generate data
      family <- test_descr$family[ii]
      model <- c("3" = "dclayton", "4" = "dgumbel", "5" = "dfrank")
      model <- model[as.character(family)]
      args <- data_sim(family = family)
      # in R
      ll_r <- VineCopula::BiCopPDF(u1 = args$udata[,1], u2 = args$udata[,2],
                                   family = family, par = args$epar, par2 = args$epar2)
      ll_r <- -sum(args$wgt * log(ll_r))
      # in TMB
      cop_adf <- TMB::MakeADFun(
        data = list(
          model = model,
          u1 = args$udata[,1],
          u2 = args$udata[,2],
          weights = args$wgt
        ),
        parameters = list(theta = args$epar),
        silent = TRUE, DLL = "LocalCop_TMBExports")
      ll_tmb <- cop_adf$fn(args$epar)
      expect_equal(ll_r, ll_tmb)
    }
  }
})



##############
# CDF TEST
##############

test_that("Copula cdf is same in VineCopula and TMB", {
  nreps <- 20
  test_descr <- expand.grid(family = c(1, 3, 4, 5), # add copula families
                            stringsAsFactors = FALSE)
  n_test <- nrow(test_descr)
  for(ii in 1:n_test) {
    for(jj in 1:nreps) {
      # generate data
      family <- test_descr$family[ii]
      model <- c("1" = "pgaussian", "3" = "pclayton", "4" = "pgumbel", "5" = "pfrank")
      model <- model[as.character(family)]
      args <- data_sim(family = family)
      # in R
      ll_r <- VineCopula::BiCopCDF(u1 = args$udata[,1], u2 = args$udata[,2],
                                   family = family, par = args$epar, par2 = args$epar2)
      ll_r <- -sum(args$wgt * log(ll_r))
      # in TMB
      cop_adf <- TMB::MakeADFun(
        data = list(
          model = model,
          u1 = args$udata[,1],
          u2 = args$udata[,2],
          weights = args$wgt
        ),
        parameters = list(theta = args$epar),
        silent = TRUE, DLL = "LocalCop_TMBExports")
      ll_tmb <- cop_adf$fn(args$epar)
      expect_equal(ll_r, ll_tmb)
    }
  }
})



############################
# PARTIAL TEST
############################


test_that("Copula partial derivative is same in VineCopula and TMB", {
  nreps <- 20
  test_descr <- expand.grid(family = c(1, 3, 4, 5), # add copula families
                            stringsAsFactors = FALSE)
  n_test <- nrow(test_descr)
  for(ii in 1:n_test) {
    for(jj in 1:nreps) {
      # generate data
      family <- test_descr$family[ii]
      model <- c("1" = "hgaussian", "3" = "hclayton", "4" = "hgumbel", "5" = "hfrank")
      model <- model[as.character(family)]
      args <- data_sim(family = family)
      # in R - VineCopula
      ll_r <- VineCopula::BiCopHfunc1(u1 = args$udata[,1], u2 = args$udata[,2], 
                                      family = family, par = args$epar, par2 = args$epar2)
      ll_r <- log(ll_r)
      ind <- ll_r > -20  # control the extremely small values in the log scale.
      ll_r <- -sum(args$wgt[ind] * ll_r[ind])
      # in TMB
      cop_adf <- TMB::MakeADFun(
        data = list(
          model = model,
          u1 = args$udata[ind,1],
          u2 = args$udata[ind,2],
          weights = args$wgt[ind]
        ),
        parameters = list(theta = args$epar[ind]),
        silent = TRUE, DLL = "LocalCop_TMBExports")
      ll_tmb <- cop_adf$fn(args$epar[ind])
      expect_equal(ll_r, ll_tmb)
      stopifnot(all.equal(ll_r, ll_tmb))
    }
  }
})

