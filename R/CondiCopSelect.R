#' Local likelihood bandwidth and/or family selection.
#'
#' Selects among a set of bandwidths and/or copula families the one which maximizes the cross-validated local likelihood.  See [CondiCopLikCV()] for details.
#'
#' @template param-u1
#' @template param-u2
#' @param family Vector of integers specifying the family set.  See [ConvertPar()].
#' @template param-x
#' @param xind Specification of `xind` for each bandwidth.  Can be a scalar integer, a vector of `nband` integers, or a list of `nband` vectors of integers.
#' @template param-degree
#' @param nu Optional vector of fixed `nu` parameter for each family.  If missing or `NA` get estimated from the data (if required)
#' @param kernel,optim_fun,cl See [CondiCopLocFit()].
#' @template param-cv_all
#' @param band Vector of positive numbers specifying the bandwidth value set.
#' @param nband If `band` is missing, automatically choose `nband` bandwidth values spanning the range of `x`.
#' @param full_out Logical; whether or not to output all fitted models or just the selected family/bandwidth combination.  See **Value**.
#' @return If `full_out = FALSE`, a list with elements `family` and `bandwidth` containing the selected value of each.  Otherwise, a list with the following elements:
#' \describe{
#'   \item{`cv`}{A data frame with `nBF = length(band) x length(family)` rows and columns named `family`, `band`, and `cv` containing the cross-validated likelihood evaluated at each combination of bandwidth and family values.}
#'   \item{`x`}{The sorted values of `x`.}
#'   \item{`eta`}{A `length(x) x nBF` matrix of eta estimates, the columns of which are in the same order as the rows of `cv`.}
#'   \item{`nu`}{A vector of length `nBF` second copula parameters, with zero if they don't exist.}
#' }
#' @example examples/CondiCopSelect.R
#' @export
CondiCopSelect <- function(u1, u2, family, x, xind = 100,
                           degree = 1, nu,
                           kernel = KernEpa, band, nband = 6,
                           optim_fun, cv_all = FALSE,
                           full_out = TRUE, cl = NA) {
  # family set
  if(missing(family)) {
    family <- .get_family(u1, u2, nper = 10)
  }
  sapply(family, .check_family)
  nfam <- length(family)
  .check_degree(degree)
  # initial parameters
  if(missing(nu)) nu <- rep(NA, nfam)
  nu <- sapply(1:nfam, function(ii) {
    .get_etaNu(u1 = u1, u2 = u2, family = family[ii],
               degree = degree, eta = c(1,0), nu = nu[ii])$nu
  })
  # bandwidth set
  if(missing(band)) band <- .get_band(x, nband)
  nband <- length(band)
  # selection process
  ## if(nband == 1 & length(family)==1) {
  ##   # no need to do selection if there is only one choice.
  ##   res <- list(band=band, family=family)
  ## } else {
  # calculate cvLIK for family & bandwidth combinations
  gridVal <- expand.grid(band = band, family = family)
  gridVal <- cbind(gridVal, nu = rep(nu, each = nband))
  # xind value adjustment for CV calculation
  if(is.numeric(xind)) {
    if(length(xind) == 1) xind <- rep(xind, nband)
    xind <- as.list(xind)
  }
  xind <- rep(xind, nfam)
  if(length(xind) != nrow(gridVal)) {
    stop("Incorrect specification of xind.")
  }
  # optimization function
  if(missing(optim_fun)) {
    optim_fun <- .optim_default
  }
  fun <- function(ii) {
    CondiCopLikCV(u1=u1, u2=u2, family = gridVal$family[ii],
                  x=x, xind = xind[[ii]], degree = degree,
                  eta=c(1,0), nu=gridVal$nu[ii], kernel=kernel,
                  band = gridVal$band[ii], optim_fun = optim_fun,
                  cveta_out = full_out, cv_all = cv_all, cl = NA)
  }
  if(!.check_parallel(cl)) {
    # run serially
    cvLIK <- sapply(1:nrow(gridVal), fun)
  } else {
    # run in parallel
    parallel::clusterExport(
      cl = cl,
      varlist = c("fun", "u1", "u2", "family", "x",
                  "band", "kernel", "optim_fun",
                  "gridVal", "xind", "cv_all",
                  "full_out"),
      envir = environment()
    )
    cvLIK <- parallel::parSapply(cl,
                                 X = 1:nrow(gridVal),
                                 FUN = fun)
  }
  if(!full_out) {
    isel <- which.max(cvLIK)
    res <- list(family = gridVal$family[isel],
                band = gridVal$band[isel])
  } else {
    res <- list(cv = cbind(gridVal[c("band", "family")],
                           cv = unlist(cvLIK["loglik",])),
                x = cvLIK["x",1]$x,
                eta = do.call(cbind, cvLIK["eta",]),
                nu = gridVal$nu)
  }
  return(res)
}

