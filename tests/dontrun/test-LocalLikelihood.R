# test C++ Local Likelihood calculations

#--- gaussian copula -----------------------------------------------------------

# R code mimics C++
GaussLL <- function(u1, u2, rho) {
  y1 <- qnorm(u1)
  y2 <- qnorm(u2)
  ll <- 2 * rho * y1 * y2
  # square everything
  y1 <- y1 * y1
  y2 <- y2 * y2
  rho <- 1 - rho*rho
  ll <- (y1 + y2 - ll)/rho + log(rho) - (y1 + y2)
  -.5 * ll
}

n <- 10
u1 <- runif(n)
u2 <- runif(n)
eta <- rnorm(1)
rho <- (exp(2*eta) - 1)/(exp(2*eta) + 1)
BiCopEta2Par(1, eta)

BiCopPDF(u1, u2, family = 1, par = rho) - exp(GaussLL(u1, u2, rho))

#--- t copula ------------------------------------------------------------------

require(VineCopula)

studentll <- function(u1, u2, rho, nu) {
  y1 <- qt(u1, nu)
  y2 <- qt(u2, nu)
  ## ll <- -(nu+2)/2 * log(1 + (y1^2 + y2^2 - 2*rho*y1*y2)/(nu*(1-rho^2)))
  ## ll <- ll + (nu+1)/2 * log(1 + y1^2/nu)
  ## ll <- ll + (nu+1)/2 * log(1 + y2^2/nu)
  ## ll <- ll - .5 * log(1-rho^2)
  ## ll
  ll <- 2 * rho * y1 * y2
  y1 <-  y1*y1
  y2 <- y2*y2
  rho <- 1 - rho*rho
  ll <- (nu + 2) * log(1.0 + (y1 + y2 - ll) / (nu * rho)) + log(rho)
  ll <- ll - ((nu + 1) * log((1.0 + y1/nu) * (1.0 + y2/nu)))
  -.5 * ll
  ##StableGammaDivision((*nu+2.0)/2.0,*nu/2.0)
  ## pow <- function(x,y) x^y
  ## ll <- 1/(nu*pi*sqrt(1.0-pow(rho,2.0))*dt(y1,nu,0)*dt(y2,nu,0))*pow(1.0+(pow(y1,2.0)+pow(y2,2.0)-2.0*rho*y1*y2)/(nu*(1.0-pow(rho,2.0))),-(nu+2.0)/2.0)
  ## log(ll)
}

n <- 10
u1 <- runif(n)
u2 <- runif(n)
eta <- rnorm(n)
rho <- (exp(2*eta) - 1)/(exp(2*eta) + 1)
nu <- rexp(1) + 5

Rcpp::sourceCpp("StudentLL.cpp")

studentll(u1, u2, rho, nu)
StudentLL(u1, u2, rho, nu)

log(BiCopPDF(u1, u2, family = 2, par = rho, par2 = nu)) - studentll(u1, u2, rho, nu)
