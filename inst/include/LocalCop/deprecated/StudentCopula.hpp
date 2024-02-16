/// @file StudentCopula.h

/// Local likelihood for the Student-t copula.
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
