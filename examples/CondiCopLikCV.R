# simulate data
family <- 5 # Frank copula
n <- 1000
X <- runif(n) # covariate values
eta_true <- 2*cos(12*pi*X) # copula dependence parameter
par_true <- BiCopEta2Par(family, eta = eta_true)
udata <- VineCopula::BiCopSim(n, family=family,
                              par = par_true$par)

# cross-validation
band <- .1
ieta <- c(1,1)
inu <- 2.5
by <- 10
cv <- CondiCopLikCV(u1= udata[,1], u2= udata[,2],
                    family=family, X=X, eta = ieta, nu = inu,
                    xind = c(seq(1, n, by = by), n), band = band,
                    cveta_out = TRUE)
