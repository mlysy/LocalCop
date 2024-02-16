/// @file GaussCopula.h

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
