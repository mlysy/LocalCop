# LocalCop: Local Likelihood Inference for Conditional Copula Models

*Elif Fidan Acar, Martin Lysy*

*March 8, 2019*

---

### Description

Implements a local likelihood estimator for the dependence parameter in bivariate conditional copula models.  Copula family and local likelihood bandwidth parameters are selected by leave-one-out cross-validation.  The models are implemented in 'TMB', meaning that the local score function is efficiently calculated via automated differentiation (AD), such that quasi-Newton algorithms may be used for parameter estimation.

### Installation

Install the **R** [**`devtools`**](https://CRAN.R-project.org/package=devtools) package and run
```{r}
devtools::install_github("mlysy/LocalCop")
```

### Examples

After installing and loading the package, please see `?CondiCopSelect`.
