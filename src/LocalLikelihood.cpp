/// @file LocalLikelihood.cpp
///
/// Local Likelihood functions.

#define TMB_LIB_INIT R_init_tsVine
#include <TMB.hpp>

/// Local likelihood for the Frank copula.
///
/// @param[in] u1 First uniform response vector.
/// @param[in] u2 Second uniform response vector.
/// @param[in] xc Covariate vector in centered form `xc = X-x`.
/// @param[in] beta Length-two vector of dependence parameters, such that `eta(xc) = beta[0] + beta[1] * xc`.
/// @param[in] wgt Local likelihood weights.
/// @return Negative local loglikelihood (scalar). 
template<class Type>
Type FrankNLL(vector<Type> u1, vector<Type> u2,
	      vector<Type> xc, vector<Type> beta, vector<Type> wgt) {
  vector<Type> theta = beta[0] + beta[1] * xc; // parameter on regular scale
  vector<Type> et1 = 1.0 - (-theta).exp(); // 1 - exp(-theta)
  vector<Type> num = theta * et1;
  vector<Type> den = et1 - (1.0 - (-theta * u1).exp()) * (1.0 - (-theta * u2).exp());
  vector<Type> ll = num / (den * den) ; // log-likelihood
  ll = wgt * (ll.log() - theta * (u1 + u2));
  return -ll.sum();
}

/// Local likelihood for the Gaussian copula.
///
/// @param[in] z1 First normal response vector.
/// @param[in] z2 Second normal response vector.
/// @param[in] xc Covariate vector in centered form `xc = X-x`.
/// @param[in] beta Length-two vector of dependence parameters, such that `eta(xc) = beta[0] + beta[1] * xc`.
/// @param[in] wgt Local likelihood weights.
/// @return Negative local loglikelihood (scalar). 
template<class Type>
Type GaussNLL(vector<Type> z1, vector<Type> z2,
	      vector<Type> xc, vector<Type> beta, vector<Type> wgt) {
  // parameter on regular scale
  // using transformation rho = (exp(eta) - exp(-eta)) / (exp(eta) + exp(-eta))
  vector<Type> rho = beta[0] + beta[1] * xc;
  rho = (2.0 * rho).exp(); 
  rho = (rho - 1.0) / (rho + 1.0);
  // start loglikelihood calculation
  vector<Type> ll = 2.0 * rho * z1 * z2;
  // precompute squares
  z1 *= z1;
  z2 *= z2;
  rho = 1.0 - rho*rho;
  ll = (z1 + z2 - ll)/rho + rho.log() - (z1 + z2);
  ll *= wgt;
  return .5 * ll.sum();
}

/// Local likelihood for the Gaussian copula.
///
/// @param[in] y1 First t-scale response vector.
/// @param[in] y2 Second t-scale response vector.
/// @param[in] xc Covariate vector in centered form `xc = X-x`.
/// @param[in] beta Length-two vector of dependence parameters, such that `eta(xc) = beta[0] + beta[1] * xc`.
/// @param[in] nu Degrees of freedom parameter (scalar).
/// @param[in] wgt Local likelihood weights.
/// @return Negative local loglikelihood (scalar). 
template<class Type>
Type StudentNLL(vector<Type> y1, vector<Type> y2,
		vector<Type> xc, vector<Type> beta, Type nu, vector<Type> wgt) {
  // parameter on regular scale
  // using transformation rho = (exp(eta) - exp(-eta)) / (exp(eta) + exp(-eta))
  vector<Type> rho = beta[0] + beta[1] * xc;
  rho = (2.0 * rho).exp(); 
  rho = (rho - 1.0) / (rho + 1.0);
  vector<Type> ll = 2.0 * rho * y1 * y2;
  // precompute squares
  y1 *= y1;
  y2 *= y2;
  rho = 1.0 - rho*rho;
  ll = (nu + 2.0) * (1.0 + (y1 + y2 - ll) / (nu * rho)).log() + rho.log();
  ll -= (nu + 1.0) * ((1.0 + y1/nu) * (1.0 + y2/nu)).log();
  ll *= wgt;
  return .5 * ll.sum();
}

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y1); // first response vector
  DATA_VECTOR(y2); // second response vector
  DATA_VECTOR(wgt); // weights
  DATA_VECTOR(xc); // centered covariates, i.e., X - x
  DATA_INTEGER(family); // copula family: 1, 2, or 5.
  PARAMETER_VECTOR(beta); // dependence parameter: eta = beta[0] + beta[1] * xc
  DATA_VECTOR(nu); // other parameter for family 2.
  Type nll = 0.0;
  if(family == 5) {
    // Frank copula
    nll = FrankNLL(y1, y2, xc, beta, wgt);
  } else if(family == 1) {
    // Gaussian copula
    nll = GaussNLL(y1, y2, xc, beta, wgt);
  } else if(family == 2) {
    // Student-t copula
    nll = StudentNLL(y1, y2, xc, beta, nu[0], wgt);
  }
  return nll;
}
