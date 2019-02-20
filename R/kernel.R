####################################
### KERNEL FUNCTIONS
####################################

# This script contains various kernel functions and the weighting scheme used in local likelihood estimation.

#' Local likelihood kernel functions.
#'
#' @name kernel
#' @aliases KernEpa KernGaus KernBeta KernBiQuad KernTriAng
#' @param t Vector of distances from mode of kernel.
#' @param par ...
#' @return Vector of kernel weights.
#' @details Describe kernels here.

#################################
# Epanechnikov kernel
#################################

#' @rdname kernel
#' @export
KernEpa <-function(t){
  return(pmax(3/4 * (1-t^2), 0))
}


#################################
# Gaussian kernel
#################################

#' @rdname kernel
#' @export
KernGaus <- function(t){
  return(dnorm(t))
}


#################################
# Beta kernel
#################################

#' @rdname kernel
#' @export
KernBeta <- function(t,par){
  return(pmax((1-t^2)^par/beta(.5,par+1), 0))
}



#################################
# BiQuadratic kernel
#################################

#' @rdname kernel
#' @export
KernBiQuad <- function(t){
  return(pmax(15/16 * (1-t^2)^2, 0))
}


#################################
# Triangular kernel
#################################

#' @rdname kernel
#' @export
KernTriAng <- function(t){
  return(pmax(1-abs(t), 0))
}
