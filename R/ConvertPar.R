###########################################
### RE-PARAMETRIZATIONS FOR COPULAS
###########################################

# This script defines functions that re-parametrize dependence parameters in bivariate copulas.
# Three types of dependence parameters are:
# 1) Calibration Function: Eta in (-Inf, Inf)
# 2) Copula Parameter: Par in family-specific range -- see BiCopParInt(family).
# 3) Kendall's Tau: Tau in [-1, 1]

#' Conversions between various bivariate copula parametrizations.
#'
#' @name ConvertPar
#' @aliases BiCopEta2Par BiCopPar2Eta BiCopEta2Tau BiCopTau2Eta
#' @param family An integer defining the bivariate copula family to use.  See **Details**.
#' @param eta,eta2 Vector of parameters on the `eta` scale.  See **Details**.
#' @param par,par2 Vector of parameters on the `par` scale.
#' @param tau Vector of parameters on the `tau` scale.
#' @details The copula family integer codes are identical to those of the \pkg{VineCopula} package. Currently, the following families are implemented:
#' \describe{
#'   \item{1}{Gaussian copula.}
#'   \item{2}{Student-t copula.}
#'   \item{3}{Clayton copula.}
#'   \item{4}{Gumbel copula.}
#'   \item{5}{Frank copula.}
#' }
#' @return Vector of converted parameters.


####################################
# Convert Eta to Par
####################################

#' @rdname ConvertPar
#' @export
BiCopEta2Par <- function(family, eta, eta2=0) {
  if(family==1 || family==2) {
    par <- (exp(eta)-exp(-eta))/(exp(eta)+exp(-eta)) # [-1,1]
    ## par2 <- eta2
  } else if(family==3 || family==13) {
    par <- exp(eta) # [0,100]
    ## par2 <- eta2
  } else if(family==4 || family==14 || family== 6 || family== 16) {
    par <- exp(eta)+1 # [1,100]
    ## par2 <- eta2
  } else if(family==5) {
    par <- eta # [-100,100]
    ## par2 <- eta2
  } else if(family==7 || family==17) {
    par <- 7*(pnorm(eta)) # [0,7]
    ## par2 <- eta2
  } else if(family==8 || family==9 || family==18 || family==19) {
    par <- 1+5*exp(eta)/(exp(eta)+1) # [1,6]
    ## par2 <- eta2
  } else if(family==10 || family==20) {
    par <- 1+7*exp(eta)/(exp(eta)+1) # [1,8]
    ## par2 <- eta2
  } else if(family==23 || family==33) {
    par <- -exp(eta) # [-100,0]
    ## par2 <- eta2
  } else if(family==24 || family==34  || family==26 || family==36) {
    par <-  -(exp(eta)+1) # [-100,-1]
    ## par2 <- eta2
  } else if(family==27 || family==37) {
    par <-  -(7*exp(eta)/(exp(eta)+1)) # [-7,0]
    ## par2 <- eta2
  } else if(family==28 || family==29 || family==38 || family==39) {
    par <-  -( 1+5*exp(eta)/(exp(eta)+1)) # [-6,-1]
    ## par2 <- eta2
  } else if(family==30 || family==40) {
    par <-  -(1+7*exp(eta)/(exp(eta)+1)) # [-8,-1]
    ## par2 <- eta2
  } else if(family==104 || family==114 || family==204 || family==214) {
    par <- exp(eta)+1 # [1,infty]
    ## par2 <- eta2
  } else if(family==124 || family==134 || family==224 || family==234) {
    par <- -(exp(eta)+1) # [-infty,-1]
    ## par2 <- eta2
  } else {
    stop("Unknown copula family.")
  }
  par2 <- eta2
  return(list(par=par, par2=par2))
}


####################################
# Convert Par to Eta
####################################

#' @rdname ConvertPar
#' @export
BiCopPar2Eta <- function(family, par, par2=0) {
  if(family==1 || family==2) {
    eta <- 0.5*log((1+par)/(1-par))
    ## eta2 <-par2
  } else if(family==3 || family==13) {
    eta <- log(par)
    ## eta2 <-par2
  } else if(family==4 || family==14 || family== 6) {
    eta <- log(par-1)
    ## eta2 <-par2
  } else if(family==5) {
    eta <- par
    ## eta2 <-par2
  } else if(family==7 || family==17) {
    eta <- qnorm(par/7)
    ## eta2 <-par2
  } else if(family==8 || family==9 || family==18 || family==19) {
    eta <- log((par-1) / (6-par))
    ## eta2 <-par2
  } else if(family==10 || family==20) {
    eta <- log((par-1) / (8-par))
    ## eta2 <-par2
  } else if(family==23 || family==33) {
    eta <- log(-par)
    ## eta2 <-par2
  } else if(family==24 || family==34  || family==26 || family==36) {
    eta <- log(-(par-1))
    ## eta2 <-par2
  } else if(family==27 || family==37) {
    eta <- log(-par/(7+par))
    ## eta2 <-par2
  } else if(family==28 || family==29 || family==38 || family==39) {
    eta <- log((-par-1) / (6+par))
    ## eta2 <-par2
  } else if(family==30 || family==40) {
    eta <- log((-par-1) / (8+par))
    ## eta2 <-par2
  } else if(family==104 || family==114 || family==204 || family==214) {
    eta <- log(par-1)
    ## eta2 <-par2
  } else if(family==124 || family==134 || family==224 || family==234) {
    eta <- log(-(par-1))
    ## eta2 <-par2
  } else {
    stop("Unknown copula family.")
  }
  eta2 <- par2
  return(list(eta = eta, eta2 = eta2))
}



####################################
# Convert Eta to Tau
####################################

#' @rdname ConvertPar
#' @export
BiCopEta2Tau <- function(family, eta, eta2=0) {
  pars <- BiCopEta2Par(family, eta, eta2)
  tau  <- VineCopula::BiCopPar2Tau(family, pars$par, pars$par2)
  return(tau)
}



####################################
# Convert Tau to Eta
####################################

#' @rdname ConvertPar
#' @export
BiCopTau2Eta <- function(family, tau) {
  par <- VineCopula::BiCopTau2Par(family,tau)
  etas <- BiCopPar2Eta(family, par, par2 = 0)
  return(etas$eta)
}

.FamilySet <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                13, 14, 16, 17, 18, 19, 20, 23,
                24, 26, 27, 28, 29, 30, 33, 34, 36,
                37, 38, 39, 40, 104, 114, 124,
                134, 204, 214, 224, 234)
