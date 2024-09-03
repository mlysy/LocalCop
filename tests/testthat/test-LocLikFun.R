#--- test local likelihood implementation in TMB -------------------------------

## library(LocalCop)
## library(TMB)
## library(testthat)
## source("helper.R")

context("LocLikFun")

test_that("LocLikFun is same in VineCopula and TMB", {
  nreps <- 20
  test_descr <- expand.grid(
    family = c(1:5, 13:14, 23:24, 33:34), # copula families
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_descr)
  for(ii in 1:n_test) {
    for(jj in 1:nreps) {
      # generate data
      family <- test_descr$family[ii]
      args <- data_sim(family = family)
      # loglik in R
      ll_r <- VineCopula::BiCopPDF(
        u1 = args$udata[,1],
        u2 = args$udata[,2],
        family = family,
        par = args$epar,
        par2 = args$epar2
      )
      ll_r <- sum(args$wgt * log(ll_r))
      # loglik in TMB
      ll_tmb <- CondiCopLocFun(
        u1 = args$udata[,1],
        u2 = args$udata[,2],
        family = family,
        x = args$x,
        x0 = args$x0,
        wgt = args$wgt,
        eta = args$eta,
        nu = args$epar2
      )
      ll_tmb <- -ll_tmb$fn(args$eta)
      expect_equal(ll_r, ll_tmb)
    }
  }
})
