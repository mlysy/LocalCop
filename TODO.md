# 29/01/2024

- [ ] Finish documentation for d/p/h of {gumbel/frank/clayton}.hpp follwing template.

- [ ] Put all of these directly on the log scale, i.e., `if(giv_log) return ans else return exp(ans)`.

- [ ] Write d/p/h for Gaussian copula in the style of the previous copulas.  This also involves checking that R and C++ version agree.  See `test/testthat` for how it was done on previous copulas.

- [ ] Write `pt()` and `qt()` in preparation for student t copula.  Note: this will be done by Martin.

- [ ] Write d/p/h for Student-t copula in style of the others.

# 27/04/2021

## Compile TMB code

After adding new TMB files, compile as follows:

```r
# setwd to within your source directory
TMBtools::export_models() # tell package to recognize TMB code has changed.
pkgbuild::compile_dll() # compile package TMB/C++ code.

# reinstall package
devtools::document() # sets up documentation and namespace
devtools::install()
```

**Tip:** Always quit + restart R at this juncture.  **TMB** C++ code is compiled into a separate `so` from the rest of the package i.e., that which **Rcpp** code lives in.  This non-standard compiling mechanism is not guaranteed to work without a fresh R session.

## R API


```r
require(LocalCop)

clayton_nll <- TMB::MakeADFun(
  # fixed inputs
  data = list(
    model = "ClaytonNLL",
    u1 = u1, 
    u2 = u2,
    weights = weights,
  ),
  # variable inputs
  parameters = list(
    # initial value to optimizer
    theta = theta0
  ),
  # config
  silent = TRUE, # for debugging
  DLL = "LocalCop_TMBExports"
)

# ClaytonNLL function and its derivatives
clayton_nll$fn(theta) # negative of weighted loglikelihood;
clayton_nll$gr(theta) # its gradient,
clayton_nll$he(theta) # hessian.
```

## Testing

The `tests/testthat` folder contains a test that can be modified to test other copula distributions. 
