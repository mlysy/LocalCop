# simulate data
family <- 5 # Frank copula
n <- 1000
X <- runif(n) # covariate values
eta_true <- 2*cos(12*pi*X) # copula dependence parameter
par_true <- BiCopEta2Par(family, eta = eta_true)
udata <- VineCopula::BiCopSim(n, family=family,
                              par = par_true$par)

# local likelihood estimation
x <- runif(20)
system.time({
  eta_hat <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                            family = family, X = X, x = x, band = .1)
})

# custom optimization routine using stats::nlminb
# (somewhat faster but less stable/accurate than default)
optim_fun <- function(obj) {
  opt <- stats::nlminb(start = obj$par, objective = obj$fn,
                       gradient = obj$gr)
  return(opt$par[1]) # always return constant term even if degree == 1
}
system.time({
  eta_hat2 <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                             family = family, X = X, x = x, band = .1,
                             optim_fun = optim_fun)
})

# relative difference
range(abs((eta_hat$eta - eta_hat2$eta)/eta_hat$eta))
