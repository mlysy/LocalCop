# LocalCop: Local Likelihood Inference for Conditional Copula Models

*Elif Fidan Acar, Martin Lysy*

*March 18, 2024*

---

### Description

Implements a local likelihood estimator for the dependence parameter in bivariate conditional copula models.  Copula family and local likelihood bandwidth parameters are selected by leave-one-out cross-validation.  The models are implemented in 'TMB', meaning that the local score function is efficiently calculated via automated differentiation (AD), such that quasi-Newton algorithms may be used for parameter estimation.

### Installation

To install the CRAN version (0.0.1):

```r
install.packages("LocalCop")
```

To install the latest development version: first install the [**`devtools`**](https://CRAN.R-project.org/package=devtools) package, then:
```{r}
devtools::install_github("mlysy/LocalCop")
```

### Usage

Please see package vignette: `vignette("LocalCop-vignette")`.
