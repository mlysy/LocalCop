#--- test parameter conversions ------------------------------------------------

context("ConvertPar")

test_that("Conversion eta => tau/par => eta gives original value.", {
  nreps <- 100
  for(ii in 1:nreps) {
    # parameters on the eta scale
    n <- sample(1:10, 1)
    etas <- list(eta = rnorm(n, sd = 1/5),
                 eta2 = runif(n, 2, 20))
    family <- sample(1:5, 1) # TODO: implement more families
    type <- sample(c("par", "tau"), 1)
    if(type == "par") {
      pars <- BiCopEta2Par(family, eta = etas$eta, eta2 = etas$eta2)
      tau <- VineCopula::BiCopPar2Tau(family, par = pars$par, par2 = pars$par2)
    } else if(type == "tau") {
      tau <- BiCopEta2Tau(family, eta = etas$eta, eta2 = etas$eta2)
      pars <- list(par = VineCopula::BiCopTau2Par(family, tau = tau),
                   par2 = etas$eta2)
    }
    # check conversions
    expect_equal(etas, BiCopPar2Eta(family, par = pars$par, par2 = pars$par2))
    expect_equal(etas$eta, BiCopTau2Eta(family, tau = tau))
  }
})
