#--- tests for student distribution --------------------------------------------

## require(TMB)
## require(numDeriv)
## require(testthat)

## model <- "pt" # this one works
## model <- "test_student" # testing
## compile(paste0(model, ".cpp"), PKG_CXXFLAGS = "-std=gnu++11")
## dyn.load(dynlib(model))

#--- pt ------------------------------------------------------------------------

# TODO: change test so that it systematically checks all cases.
test_that("pt gives same value in R and TMB", {
  ntest <- 100
  sim_arg <- function(n) {
    x <- sample(c(1, 10, 100), n, replace = TRUE) + rnorm(n)
    sample(c(-1, 1), n, replace = TRUE) * x
  }
  for(ii in 1:ntest) {
    n <- sample(c(1, 2, 10), 1)
    pt_adf <- TMB::MakeADFun(data = list(model = "pt"),
                             parameters = list(q = rep(0, n), df = rep(1, n)),
                             silent = TRUE,
                             DLL = "LocalCop_TMBExports")
    q <- sim_arg(n)
    df <- abs(sim_arg(n))
    p_tmb <- pt_adf$fn(c(q, df))
    ## suppressWarnings({
    p_r <- sum(pt(q, df))
    ## })
    expect_equal(p_r, p_tmb)
  }
})

#--- qt ------------------------------------------------------------------------


test_that("qt gives same value in R and TMB", {
  ntest <- 100
  sim_arg <- function(n) {
    x <- sample(c(1, 10, 100), n, replace = TRUE) + rnorm(n)
    sample(c(-1, 1), n, replace = TRUE) * x
  }
  error <- function(x1, x2) {
    abs(x1 - x2)/(.5 * abs(x1 + x2) + .1)
  }
  for(ii in 1:ntest) {
    n <- sample(c(1, 2, 10), 1)
    qt_adf <- TMB::MakeADFun(data = list(model = "qt"),
                             parameters = list(p = rep(.5, n), df = rep(1, n)),
                             silent = TRUE)
    p <- runif(n)
    df <- abs(sim_arg(n))
    q_tmb <- qt_adf$fn(c(p, df))
    ## suppressWarnings({
    q_r <- sum(qt(p, df))
    ## })
    expect_equal(error(q_r, q_tmb) < 1e-10, TRUE)
  }
})
