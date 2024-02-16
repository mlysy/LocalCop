/// @file LocalLikelihood.hpp
///
/// Local Likelihood functions.

#include "LocalCop/deprecated/FrankCopula.hpp"
#include "LocalCop/deprecated/GaussCopula.hpp"
#include "LocalCop/deprecated/GumbelCopula.hpp"
#include "LocalCop/deprecated/StudentCopula.hpp"
#include "LocalCop/deprecated/ClaytonCopula.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type LocalLikelihood(objective_function<Type> *obj) {
  DATA_VECTOR(y1); // first response vector
  DATA_VECTOR(y2); // second response vector
  DATA_VECTOR(wgt); // weights
  DATA_VECTOR(xc); // centered covariates, i.e., X - x
  DATA_INTEGER(family); // copula family: 1-5.
  PARAMETER_VECTOR(beta); // dependence parameter: eta = beta[0] + beta[1] * xc
  DATA_VECTOR(nu); // other parameter for family 2.
  Type nll = 0.0;
  if(family == 1) {
    // Gaussian copula
    nll = GaussNLL(y1, y2, xc, beta, wgt);
  } else if(family == 2) {
    // Student-t copula
    nll = StudentNLL(y1, y2, xc, beta, nu[0], wgt);
  } else if(family == 3) {
    // Clayton copula
    nll = ClaytonNLL(y1, y2, xc, beta, wgt);
  } else if(family == 4) {
    // Gumbel copula
    nll = GumbelNLL(y1, y2, xc, beta, wgt);
  } else if(family == 5) {
    // Frank copula
    nll = FrankNLL(y1, y2, xc, beta, wgt);
  }
  return nll;
}
