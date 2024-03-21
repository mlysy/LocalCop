# LocalCop 0.0.1

- Initial CRAN submission.

## Major Changes

- **Breaking change:** replaced function arguments `X` with `x` and `x` with `x0`.

- Using **TMBTools** method to export **TMB** models.

- Rewrote  `LocalLikelihood` **TMB** model template to use copula densities `dGaussian()`, `dClayton()`, etc., located in `inst/include/LocalCop`.  This involved writing `pt()` and `qt()` functions in **TMB** for Student-t copula `dStudent()`.

- Added vignette.

## Minor Changes

- Using Markdown in **roxygen** documentation.

- Refactored tests, including fixing random seed.

- Fixed seed in `CondiCopSelect()` example.

- Fixed plot in `CondiCopLocFit()` example.

# LocalCop 0.0.0.9000

- Initial GitHub public submission.
