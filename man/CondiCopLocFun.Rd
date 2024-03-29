% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CondiCopLocFun.R
\name{CondiCopLocFun}
\alias{CondiCopLocFun}
\title{Create a \pkg{TMB} local likelihood function.}
\usage{
CondiCopLocFun(u1, u2, family, x, x0, wgt, degree = 1, eta, nu)
}
\arguments{
\item{u1}{Vector of first uniform response.}

\item{u2}{Vector of second uniform response.}

\item{family}{An integer defining the bivariate copula family to use.  See \code{\link[=ConvertPar]{ConvertPar()}}.}

\item{x}{Vector of observed covariate values.}

\item{x0}{Scalar covariate value at which to evaluate the local likelihood.  Does not have to be a subset of \code{x}.}

\item{wgt}{Vector of positive kernel weights.}

\item{degree}{Integer specifying the polynomial order of the local likelihood function.  Currently only 0 and 1 are supported.}

\item{eta}{Value of the copula dependence parameter.  Scalar or vector of length two, depending on whether \code{degree} is 0 or 1.}

\item{nu}{Value of the other copula parameter.  Scalar or vector of same length as \code{u1}.  Ignored if \code{family != 2}.}
}
\value{
A list as returned by a call to \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}.  In particular, this contains elements \code{fun} and \code{gr} for the \emph{negative} local likelihood and its gradient with respect to \code{eta}.
}
\description{
Wraps a call to \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}.
}
\examples{
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
                          x = rep(0, n), x0 = 0, # centered covariate x - x0 == 0
                          wgt = rep(1, n), # unweighted
                          degree = 0, # zero-order fit
                          eta = c(rho2eta(rho), 0),
                          nu = nu)

# likelihood function: recall that TMB requires a _negative_ ll
stucop_lik <- function(rho) {
  -nll_obj$fn(c(rho2eta(rho), 0))
}

# compare to VineCopula.
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
range(ll1 - ll2)
}
