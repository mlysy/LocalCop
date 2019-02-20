#' Local likelihood estimation at a single covariate value.
#'
#' @param u1 Vector of first uniform response.
#' @param u2 Vector of second uniform response.
#' @template param-family
#' @param z Vector of centered covariate values, i.e., \code{X-x}.
#' @param wgt Vector of kernel weights.
#' @param degree Character vector "constant" or "linear" specifying the type of local likelihood fit.
#' @param eta Initial value of the copula dependence parameter.
#' @param nu Value of other copula parameters (if they exist).
#' @return List with elements \code{eta} and \code{loglik}.
#' @export
CondiCopLocFit1 <- function(u1, u2, family,
                            z, wgt, degree = c("linear", "constant"),
                            eta, nu) {
  obj <- CondiCopLocFun(u1 = u1, u2 = u2, family = family,
                        z = z, wgt = wgt, degree = degree, eta = eta, nu = nu)
  opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
               method = "BFGS")
  return(list(eta = as.numeric(opt$par), loglik = -opt$value))
}

#--- helper function (not exported) --------------------------------------------

.format_degree <- function(degree = c("linear", "constant")) {
  degree <- match.arg(degree)
  return(as.numeric(degree == "linear"))
}

#' Create a TMB local likelihood function.
#'
#' Wraps a call to \code{TMB::MakeADFun}.
#'
#' @param u1 Vector of first uniform response.
#' @param u2 Vector of second uniform response.
#' @template param-family
#' @param z Vector of centered covariate values, i.e., \code{X-x}.
#' @param wgt Vector of kernel weights.
#' @param degree Character vector "constant" or "linear" specifying the type of local likelihood fit.
#' @param eta Initial value of the copula dependence parameter.
#' @param nu Value of other copula parameters (if they exist).
#' @return A list as returned by a call to \code{TMB::MakeADFun}.
#' @export
CondiCopLocFun <- function(u1, u2, family,
                           z, wgt, degree = c("linear", "constant"),
                           eta, nu) {
  wpos <- wgt > 0 # index of positive weights
  # create TMB function: need data, parameters,
  # and map (for fixing eta[2] if needed)
  odata <- list(u1 = u1[wpos], u2 = u2[wpos],
                wgt = wgt[wpos], z = z[wpos]) # TODO: family, nu
  oparam <- list(beta = eta)
  degree <- .format_degree(degree) # convert degree to TMB::map
  omap <- list(beta = factor(c(1, 2)))
  if(degree == 0) omap$beta[2] <- NA
  TMB::MakeADFun(data = odata, parameters = oparam,
                 map = omap, DLL = "tsVine", silent = TRUE)
}
