/// @file LocalLikelihood.hpp
///
/// @brief Local Likelihood calculations for the five major families.

#include "LocalCop/frank.hpp"
#include "LocalCop/gaussian.hpp"
#include "LocalCop/gumbel.hpp"
#include "LocalCop/student.hpp"
#include "LocalCop/clayton.hpp"

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
  vector<Type> theta = beta(0) + beta(1) * xc;
  vector<Type> lpdf(theta.size());
  if(family == 1) {
    // Gaussian copula
    theta = (2.0 * theta).exp(); 
    theta = (theta - 1.0) / (theta + 1.0);
    lpdf = LocalCop::dgaussian(y1, y2, theta, 1);
  } else if(family == 2) {
    // Student-t copula
    theta = (2.0 * theta).exp(); 
    theta = (theta - 1.0) / (theta + 1.0);
    lpdf = LocalCop::dstudent(y1, y2, theta, nu, 1);
  } else if(family == 3) {
    // Clayton copula
    theta = theta.exp();
    lpdf = LocalCop::dclayton(y1, y2, theta, 1);
  } else if(family == 4) {
    // Gumbel copula
    theta = 1.0 + theta.exp();
    lpdf = LocalCop::dgumbel(y1, y2, theta, 1);
  } else if(family == 5) {
    // Frank copula
    lpdf = LocalCop::dfrank(y1, y2, theta, 1);
  } else {
    Rf_error("Unknown copula family.");
  }
  lpdf.array() *= wgt.array();
  nll = -sum(lpdf);
  return nll;
}
