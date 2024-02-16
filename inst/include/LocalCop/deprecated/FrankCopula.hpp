/// @file FrankCopula.h

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
