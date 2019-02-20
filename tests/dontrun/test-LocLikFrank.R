#--- test the Frank Local likelihood model --------------------------------

require(tsVine)
require(testthat)
# run from any subfolder of tsVine
source(devtools::package_file("tests", "testthat", "loclik-functions.R"))

# generate data
set.seed(180219)

n <- 1000 # number of time points
X <- (1:n)/n  # time points in  (0,1] interval
family <- 5   # Frank copula
eta.fnc <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
true.par <- BiCopEta2Par(family = family, eta.fnc(X))$par
true.tau <- BiCopEta2Tau(family = family, eta.fnc(X))
udata <- VineCopula::BiCopSim(N=n, family=family, par = true.par)

# function value
eta <- c(1,1)
band <- .1
x0 <- mean(X)
wgt <- KernWeight(X = X, x = x0, band = band)

# R
fn_r <- function(eta) {
  -llfrank(u1 = udata[,1],
           u2 = udata[,2],
           z = X - x0, wgt = wgt, eta = eta)
}

# TMB
obj_tmb <- CondiCopLocFun(u1 = udata[,1], u2 = udata[,2],
                          family = 5, z = X-mean(X), wgt = wgt,
                          eta = eta, nu = 0)

test_that("R and TMB local likelihoods are the same.", {
  nll_r <- fn_r(eta)
  nll_tmb <- obj_tmb$fn(eta)
  expect_equal(nll_tmb, nll_r)
})

test_that("R and TMB optimization are the same.", {
  opt_r <- optim(par = eta, fn = fn_r)
  eta_r <- as.numeric(opt_r$par)
  ll_r <- -opt_r$value
  opt_tmb <- CondiCopLocFit1(u1 = udata[,1], u2 = udata[,2],
                             family = 5, z = X-mean(X), wgt = wgt,
                             eta = eta, nu = 0)
  eta_tmb <- as.numeric(opt_tmb$eta)
  ll_tmb <- opt_tmb$loglik
  expect_equal(eta_r, eta_tmb, tol = .01,
               scale = if(ll_r > ll_tmb) abs(eta_r) else abs(eta_tmb))
  expect_equal(ll_r, ll_tmb, tol = .01,
               scale = if(ll_r > ll_tmb) abs(ll_r) else abs(ll_tmb))
})

#--- parallel check ------------------------------------------------------------

require(parallel)
cl <- makeCluster(4)

# check multiple version
system.time({
  opt1 <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                         family = 5, X = X,
                         eta = eta, nu = 0, band = .1, cl = NA)
})

system.time({
  opt2 <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                         family = 5, X = X + rnorm(1),
                         eta = eta, nu = 0, band = .1, cl = cl)
})

stopCluster(cl)
