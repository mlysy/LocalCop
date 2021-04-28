/// @file ClaytonCopula.h

/// Local likelihood for the Clayton copula.
///
/// @param[in] lu1 First log-uniform response vector.
/// @param[in] lu2 Second log-uniform response vector.
/// @param[in] xc Covariate vector in centered form `xc = X-x`.
/// @param[in] beta Length-two vector of dependence parameters, such that `eta(xc) = beta[0] + beta[1] * xc`.
/// @param[in] wgt Local likelihood weights.
/// @return Negative local loglikelihood (scalar). 
template<class Type>
Type ClaytonNLL(vector<Type> lu1, vector<Type> lu2,
	      vector<Type> xc, vector<Type> beta, vector<Type> wgt) {
  // parameter on regular scale
  vector<Type> theta = (beta[0] + beta[1] * xc).exp();
  vector<Type> ll = (1.0 + theta).log() - (1.0 + theta) * (lu1 + lu2);
  ll -= (2.0 + 1.0/theta) * ((-theta * lu1).exp() + (-theta * lu2).exp() - 1.0).log();
  ll *= wgt;
  return -ll.sum();
}
