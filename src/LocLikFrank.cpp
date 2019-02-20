// Local Likelihood for Frank copula -- vectorized version

#define TMB_LIB_INIT R_init_tsVine
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(u1); // first uniform response
  DATA_VECTOR(u2); // second uniform response
  DATA_VECTOR(wgt); // weights
  DATA_VECTOR(z); // "centered" covariates, i.e., xi - x0
  PARAMETER_VECTOR(beta); // theta = eta = beta[0] + beta[1] * z
  vector<Type> theta = beta[0] + beta[1] * z; // parameter on regular scale
  vector<Type> et1 = 1.0 - (-theta).exp(); // 1 - exp(-theta)
  vector<Type> num = theta * et1;
  vector<Type> den = et1 - (1.0 - (-theta * u1).exp()) * (1.0 - (-theta * u2).exp());
  vector<Type> ll = num / (den * den) ; // log-likelihood
  ll = wgt * (ll.log() - theta * (u1 + u2));
  return -ll.sum();
}
