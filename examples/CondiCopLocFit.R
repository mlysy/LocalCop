# simulate data
family <- 5 # Frank copula
n <- 1000
X <- runif(n) # covariate values
eta_fun <- function(x) 2*cos(12*pi*x) # copula dependence parameter
eta_true <- eta_fun(X)
par_true <- BiCopEta2Par(family, eta = eta_true)
udata <- VineCopula::BiCopSim(n, family=family,
                              par = par_true$par)

# local likelihood estimation
x <- seq(min(X), max(X), len = 100)
band <- .02
system.time({
  eta_hat <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                            family = family, X = X, x = x, band = band)
})

# custom optimization routine using stats::nlminb
# (somewhat faster but less stable/accurate than default)
optim_fast <- function(obj) {
  opt <- stats::nlminb(start = obj$par, objective = obj$fn,
                       gradient = obj$gr)
  return(opt$par[1]) # always return constant term even if degree == 1
}
system.time({
  eta_hat2 <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                             family = family, X = X, x = x, band = band,
                             optim_fun = optim_fast)
})

plot(x, BiCopEta2Tau(family, eta = eta_fun(x)), type = "l",
     xlab = expression(x), ylab = expression(tau(x)))
lines(x, BiCopEta2Tau(family, eta = eta_hat$eta), col = "red")
lines(x, BiCopEta2Tau(family, eta = eta_hat2$eta), col = "blue")
legend("bottomright", legend = c("True", "optim_default", "optim_fast"))
