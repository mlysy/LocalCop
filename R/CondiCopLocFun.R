#' Create a \pkg{TMB} local likelihood function.
#'
#' Wraps a call to \code{TMB::MakeADFun}.
#'
#' @template param-u1
#' @template param-u2
#' @template param-family
#' @template param-X
#' @param x Scalar covariate value at which to evaluate the local likelihood.  Does not have to be a subset of \code{X}.
#' @param wgt Vector of positive kernel weights.
#' @template param-degree
#' @param eta Value of the copula dependence parameter.  Scalar or vector of length two, depending on whether \code{degree} is 0 or 1.
#' @param nu Value of the other copula parameter.  Needn't be supplied if it doesn't exist.
#' @return A list as returned by a call to \code{TMB::MakeADFun}.  In particular, this contains elements \code{fun} and \code{gr} for the *negative* local likelihood and its gradient with respect to \code{eta}.
#' @export
CondiCopLocFun <- function(u1, u2, family,
                           X, x, wgt, degree = 1,
                           eta, nu) {
  if(!family %in% 1:5) {
    stop("Unsupported copula family (must be integer between 1-5).")
  }
  wpos <- wgt > 0 # index of positive weights
  # create TMB function
  # data input
  if(family != 2) nu <- 0 # second copula parameter
  odata <- list(y1 = u1[wpos], y2 = u2[wpos],
                wgt = wgt[wpos], xc = X[wpos]-x,
                family = family, nu = nu)
  # pre-transform data to correct scale
  if(family == 1) {
    # Gaussian copula
    odata$y1 <- qnorm(odata$y1)
    odata$y2 <- qnorm(odata$y2)
  }
  if(family == 2) {
    # Student-t copula
    odata$y1 <- qt(odata$y1, df = nu)
    odata$y2 <- qt(odata$y2, df = nu)
  }
  if(family == 3) {
    # Clayton copula
    odata$y1 <- log(odata$y1)
    odata$y2 <- log(odata$y2)
  }
  if(family == 4) {
    # Gumbel copula
    odata$y1 <- log(-log(odata$y1))
    odata$y2 <- log(-log(odata$y2))
  }
  # 2nd copula parameter
  ## if(family != 2) odata$nu <- 0
  oparam <- list(beta = eta)
  # convert degree to TMB::map
  ## degree <- .format_degree(degree)
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  omap <- list(beta = factor(c(1, 2)))
  if(degree == 0) {
    oparam$beta[2] <- 0
    omap$beta[2] <- NA
  }
  TMB::MakeADFun(data = odata, parameters = oparam,
                 map = omap, DLL = "LocalCop", silent = TRUE)
}

## #' Local likelihood estimation at a single covariate value.
## #'
## #' @inheritParams CondiCopLocFun
## #' @return List with elements \code{eta} and \code{loglik}.
## #' @export
## CondiCopLocFit1 <- function(u1, u2, family,
##                             X, x, wgt, degree = c("linear", "constant"),
##                             eta, nu = 0) {
##   obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
##                         z = z, wgt = wgt, degree = degree, eta = eta, nu = nu)
##   opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
##                method = "BFGS")
##   return(list(eta = as.numeric(opt$par), loglik = -opt$value))
## }

#--- helper function (not exported) --------------------------------------------

# convert degree to integer
.format_degree <- function(degree) {
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  ## degree <- match.arg(degree)
  ## return(as.numeric(degree == "linear"))
}
