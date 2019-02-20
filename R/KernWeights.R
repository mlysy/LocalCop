#################################
# Kernel Weights
#################################

#' Calculate local likelihood kernel weights.
#'
#' @template param-X
#' @template param-x0
#' @param band Bandwidth parameter (positive scalar).
#' @export
KernWeight <- function(X, x, band, kernel = KernEpa, band.method ="constant"){

  stopifnot(length(band)==1 | length(x)==1)

  if(band.method=="constant"){
    hval <- band
    # w <- outer(X, x, function(Y,y) (1/band)*kernel((Y-y)/band))
    }

  if(band.method=="variable"){
    if(band > 1) stop("band should be <= 1 for the variable bandwidth method.")
    k <- as.integer(band*length(X))
    if(band == 1){k <- k-1}
    hval <- max(sort(abs(X-x))[1:(k+1)])
   # hval <-  sapply(x, function(y) max(sort(abs(X-y))[1:(k+1)]))
   # w <- sapply(1:length(x), function(k)(1/hval[k])*kernel((X-x[k])/hval[k]) )
  }

  w <- (1/hval)*kernel((X-x)/hval)

  # w <- sapply(hval, function(h) outer(X, x, function(Y,y) (1/h)*kernel((Y-y)/h)), simplify = "array")
  return(w)
}



