# simulate data
family <- 5 # Frank copula
n <- 1000
X <- runif(n) # covariate values
eta_fun <- function(x) 2*cos(12*pi*x) # copula dependence parameter
eta_true <- eta_fun(X)
par_true <- BiCopEta2Par(family, eta = eta_true)
udata <- VineCopula::BiCopSim(n, family=family,
                              par = par_true$par)

# bandwidth and family selection
bandset <- c(.01, .02, .05, .1)
famset <- c(2, 5)
xind <- 200
system.time({
  cvsel <- CondiCopSelect(u1= udata[,1], u2 = udata[,2],
                          X = X, family = famset, band = bandset,
                          xind = xind)
})

# compare estimates to true value
xseq <- cvsel$x
famsel <- cvsel$cv$family
bandsel <- cvsel$cv$band
etasel <- cvsel$eta
clrs <- c("red", "blue", "orange", "green4")
names(clrs) <- bandset
ltys <- 1:2
names(ltys) <- famset
par(mfrow = c(1,2))
# student t copula
plot(xseq, BiCopEta2Tau(family, eta = eta_fun(xseq)),
     type = "l", lwd = 2, ylim = c(-.5, .5))
for(ii in 1:length(bandset)) {
  lines(xseq, BiCopEta2Tau(famsel[ii], eta = etasel[,ii]),
        col = clrs[as.character(bandsel[ii])],
        lty = ltys[as.character(famsel[ii])], lwd = 1)
}
# frank copula
plot(xseq, BiCopEta2Tau(family, eta = eta_fun(xseq)),
     type = "l", lwd = 2, ylim = c(-.5, .5))
for(ii in length(bandset)+1:length(bandset)) {
  lines(xseq, BiCopEta2Tau(famsel[ii], eta = etasel[,ii]),
        col = clrs[as.character(bandsel[ii])],
        lty = ltys[as.character(famsel[ii])], lwd = 1)
}
