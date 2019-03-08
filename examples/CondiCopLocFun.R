# the following example shows how to create
# an unconditional copula likelihood function

# simulate data
n <- 1000 # sample size
family <- 2 # Student-t copula
rho <- runif(1, -1, 1) # unconditional dependence parameter
nu <- runif(1, 4, 20)# degrees of freedom parameter
udata <- VineCopula::BiCopSim(n, family = family, par = rho, par2 = nu)

# create likelihood function

# parameter conversion: equivalent to BiCopPar2Eta(family = 2, ...)
rho2eta <- function(rho) .5 * log((1+rho)/(1-rho))
nll_obj <- CondiCopLocFun(u1 = udata[,1], u2 = udata[,2], family = family,
                          X = rep(0, n), x = 0, # centered covariate X - x == 0
                          wgt = rep(1, n), # unweighted
                          degree = 0, # zero-order fit
                          eta = c(rho2eta(rho), 0),
                          nu = nu)

# likelihood function: recall that TMB requires a _negative_ ll
stucop_lik <- function(rho) {
  -nll_obj$fn(c(rho2eta(rho), 0))
}

# compare to VineCopula.  since LocalCop considers nu to be fixed, its student-t
# likelihood differs from that of VineCopula in two ways:
# 1. qt(u1, df = nu) and qt(u2, df = nu) only need to be calculated once.
# 2. nu-only terms are dropped form LocalCop likelihood.
# these factors contribute to the difference in speed observed below.
rhovec <- runif(50, -1, 1)
system.time({
  ll1 <- sapply(rhovec, stucop_lik) # LocalCop
})
system.time({
  ll2 <- sapply(rhovec, function(rho) {
    # VineCopula
    sum(log(VineCopula::BiCopPDF(u1 = udata[,1], u2 = udata[,2],
                                        family = family,
                                        par = rho, par2 = nu)))
  })
})

# difference between the two
range(ll1 -
      (ll2 - n * (lgamma(.5*(nu+2)) + lgamma(.5*nu) - 2*lgamma(.5*(nu+1)))))
