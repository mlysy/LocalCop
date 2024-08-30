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
  // rotated copulas
  vector<Type> u1 = y1;
  vector<Type> u2 = y2;
  int fam = family;
  if((family == 13) | (family == 14)) {
    // 180 degree rotation
    u1 = Type(1.0) - u1;
    u2 = Type(1.0) - u2;
    fam = family - 10;
  }
  if((family == 23) | (family == 24)) {
    // 90 degree rotation
    u1 = Type(1.0) - u1;
    fam = family - 20;
  }
  if((family == 33) | (family == 34)) {
    // 270 degree rotation
    u2 = Type(1.0) - u2;
    fam = family - 30;
  }
  Type nll = 0.0;
  vector<Type> theta = beta(0) + beta(1) * xc;
  vector<Type> lpdf(theta.size());
  if(fam == 1) {
    // Gaussian copula
    theta = (2.0 * theta).exp(); 
    theta = (theta - 1.0) / (theta + 1.0);
    lpdf = LocalCop::dgaussian(u1, u2, theta, 1);
  } else if(fam == 2) {
    // Student-t copula
    theta = (2.0 * theta).exp(); 
    theta = (theta - 1.0) / (theta + 1.0);
    lpdf = LocalCop::dstudent(u1, u2, theta, nu, 1);
  } else if(fam == 3) {
    // Clayton copula
    theta = theta.exp();
    lpdf = LocalCop::dclayton(u1, u2, theta, 1);
  } else if(fam == 4) {
    // Gumbel copula
    theta = 1.0 + theta.exp();
    lpdf = LocalCop::dgumbel(u1, u2, theta, 1);
  } else if(fam == 5) {
    // Frank copula
    lpdf = LocalCop::dfrank(u1, u2, theta, 1);
  } else {
    Rf_error("Unknown copula family.");
  }
  lpdf.array() *= wgt.array();
  nll = -sum(lpdf);
  return nll;
}
