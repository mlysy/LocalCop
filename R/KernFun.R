####################################
### KERNEL FUNCTIONS
####################################

# This script contains various kernel functions and the weighting scheme used in local likelihood estimation.

#' Local likelihood kernel functions.
#'
#' @name KernFun
#' @aliases KernEpa KernGaus KernBeta KernBiQuad KernTriAng
#' @param t Vector of distances from mode of kernel.
#' @param par Shape parameter for Beta kernel (positive scalar).
#' @return Vector of kernel weights.
#' @details Describe kernels here.
NULL

#################################
# Epanechnikov kernel
#################################

#' @rdname KernFun
#' @export
KernEpa <-function(t) {
  return(pmax(3/4 * (1-t^2), 0))
}


#################################
# Gaussian kernel
#################################

#' @rdname KernFun
#' @export
KernGaus <- function(t) {
  return(dnorm(t))
}


#################################
# Beta kernel
#################################

#' @rdname KernFun
#' @export
KernBeta <- function(t,par=.5) {
  return(pmax(1-t^2, 0)^par/beta(.5,par+1))
}



#################################
# BiQuadratic kernel
#################################

#' @rdname KernFun
#' @export
KernBiQuad <- function(t) {
  return(15/16 * pmax(1-t^2, 0)^2)
}


#################################
# Triangular kernel
#################################

#' @rdname KernFun
#' @export
KernTriAng <- function(t) {
  return(pmax(1-abs(t), 0))
}
