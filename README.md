
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LocalCop: Local Likelihood Inference for Conditional Copula Models

*Elif Fidan Acar, Martin Lysy*

<!-- badges: start -->

[![R-CMD-check](https://github.com/mlysy/LocalCop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mlysy/LocalCop/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

------------------------------------------------------------------------

### Description

Implements a local likelihood estimator for the dependence parameter in
bivariate conditional copula models. Copula family and local likelihood
bandwidth parameters are selected by leave-one-out cross-validation. The
models are implemented in [**TMB**](https://github.com/kaskr/adcomp),
meaning that the local score function is efficiently calculated via
automated differentiation (AD), such that quasi-Newton algorithms may be
used for parameter estimation.

### Installation

To install the CRAN version (0.0.1):

``` r
install.packages("LocalCop", INSTALL_opts = "--install-tests")
```

To install the latest development version: first install the
[**devtools**](https://CRAN.R-project.org/package=devtools) package,
then:

``` r
devtools::install_github("mlysy/LocalCop", INSTALL_opts = "--install-tests")
```

### Usage

Please see package vignette: `vignette("LocalCop-vignette")`.

### Unit Tests

To verify that the package has been installed correctly, you can run its
unit tests. First install the
[**testthat**](https://CRAN.R-project.org/package=testthat) package,
then:

``` r
testthat::test_package("LocalCop", reporter = "progress")
```

### Contributing

Contributions in the form of bug reports, fixes, extensions,
improvements, etc. are most welcome. Please file an
[issue](https://github.com/mlysy/LocalCop/issues) before submitting a
PR.
