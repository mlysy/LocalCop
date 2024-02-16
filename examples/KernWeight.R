X <- sort(runif(20))
x <- runif(1, min = min(X), max= max(X))
KernWeight(X, x, band=0.3, kernel = KernEpa, band_type = "constant")
KernWeight(X, x, band=0.3, kernel = KernEpa, band_type = "variable")
