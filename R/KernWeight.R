#################################
# Kernel Weights
#################################

#' Calculate local likelihood kernel weights.
#'
#' @template param-X
#' @template param-x0
#' @template param-kernel
#' @param band Kernel bandwidth parameter (positive scalar).  See \strong{Details.}
#' @param band_type  A character string specifying the type of bandwidth: either "constant" or "variable".  See \strong{Details}.
#' @details For the constant bandwidth of size \code{band = h}, the weights are calculated as
#' \preformatted{
#' wgt = kernel((X-x) / h) / h
#' }
#' where \code{kernel} is the kernel function.  For bandwidth type "variable", a fixed fraction \code{band} of observations is used, i.e,
#' \preformatted{
#' h = sort( abs(X-x) )[ floor(band*length(X)) ]
#' }
#' @example
#' X <- sort(runif(20))
#' x <- runif(1, min = min(X), max= max(X))
#' KernWeight(X, x, band=0.3, kernel = KernEpa, band_type = "constant") 
#' KernWeight(X, x, band=0.3, kernel = KernEpa, band_type = "variable") 
#' 
#' @export
KernWeight <- function(X, x, band, kernel = KernEpa, band_type = "constant") {
  stopifnot(length(band)==1 | length(x)==1)
  band_type <- match.arg(band_type, choices = c("constant", "variable"))
  if(band_type=="constant") {
    hval <- band
    # w <- outer(X, x, function(Y,y) (1/band)*kernel((Y-y)/band))
  } else if(band_type=="variable") {
    if(band > 1) stop("band should be <= 1 for the variable bandwidth method.")
    k <- as.integer(band*length(X))
    if(band == 1){k <- k-1}
    hval <- max(sort(abs(X-x))[1:(k+1)])
   # hval <-  sapply(x, function(y) max(sort(abs(X-y))[1:(k+1)]))
   # w <- sapply(1:length(x), function(k)(1/hval[k])*kernel((X-x[k])/hval[k]) )
  }
  # weights
  w <- (1/hval)*kernel((X-x)/hval)
  # w <- sapply(hval, function(h) outer(X, x, function(Y,y) (1/h)*kernel((Y-y)/h)), simplify = "array")
  return(w)
}

