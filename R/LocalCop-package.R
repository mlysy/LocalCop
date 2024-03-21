#' Local likelihood inference for conditional copula models.
#'
#' Fits a bivariate conditional copula \eqn{C(u_1, u_2 | \theta_x)}, where \eqn{\theta_x} is a variable dependence parameter, nonparametrically estimated from a single covariate \eqn{x} via local likelihood.
#' @example examples/CondiCopSelect.R
#' @rawNamespace useDynLib(LocalCop, .registration=TRUE); useDynLib(LocalCop_TMBExports)
#' @importFrom stats approx dnorm optim pnorm qnorm qt cor
"_PACKAGE"
