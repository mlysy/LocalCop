require(tsVine)

# generate data
set.seed(180219)

n <- 1000 # number of time points
X <- (1:n)/n  # time points in  (0,1] interval
family <- 5   # Frank copula
eta.fnc <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
true.par <- BiCopEta2Par(family = family, eta.fnc(X))$par
true.tau <- BiCopEta2Tau(family = family, eta.fnc(X))
udata <- VineCopula::BiCopSim(N=n, family=family, par = true.par)

CondiCopLikCV(u1 = udata[,1], u2 = udata[,2], family = family,
              X = X, xind = c(seq(1, n, by = 20), n), band = .1)
