#################################
# Kernel Weights
#################################

#' Calculate local likelihood kernel weights.
#'
#' @template param-x
#' @template param-x0
#' @template param-kernel
#' @param band Kernel bandwidth parameter (positive scalar).  See **Details.**
#' @param band_type  A character string specifying the type of bandwidth: either "constant" or "variable".  See **Details**.
#' @return A vector of nonnegative kernel weights of the same length as `x`.
#' @details For the constant bandwidth of size `band = h`, the weights are calculated as
#' ```
#' wgt = kernel((x-x0) / h) / h
#' ```
#' where `kernel` is the kernel function.  For bandwidth type "variable", a fixed fraction `band` of observations is used, i.e,
#' ```
#' h = sort( abs(x-x0) )[ floor(band*length(x)) ]
#' ```
#' @example examples/KernWeight.R
#' @export
KernWeight <- function(x, x0, band, kernel = KernEpa, band_type = "constant") {
  stopifnot(length(band)==1 | length(x0)==1)
  band_type <- match.arg(band_type, choices = c("constant", "variable"))
  if(band_type=="constant") {
    hval <- band
    # w <- outer(x, x0, function(Y,y) (1/band)*kernel((Y-y)/band))
  } else if(band_type=="variable") {
    if(band > 1) stop("band should be <= 1 for the variable bandwidth method.")
    k <- as.integer(band*length(x))
    if(band == 1) k <- k-1
    hval <- max(sort(abs(x-x0))[1:(k+1)])
    # hval <-  sapply(x0, function(y) max(sort(abs(x-y))[1:(k+1)]))
    # w <- sapply(1:length(x0), function(k)(1/hval[k])*kernel((x-x0[k])/hval[k]) )
  }
  # weights
  w <- (1/hval)*kernel((x-x0)/hval)
  # w <- sapply(hval, function(h) outer(x, x0, function(Y,y) (1/h)*kernel((Y-y)/h)), simplify = "array")
  return(w)
}

