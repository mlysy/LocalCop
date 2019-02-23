/// @file GumbelCopula.h

/// Local likelihood for the Gumbel copula.
///
/// @param[in] lu1 First loglog-uniform response vector, i.e., 'lu1 = log(-log(u1))'.
/// @param[in] lu2 Second loglog-uniform response vector.
/// @param[in] xc Covariate vector in centered form `xc = X-x`.
/// @param[in] beta Length-two vector of dependence parameters, such that `eta(xc) = beta[0] + beta[1] * xc`.
/// @param[in] wgt Local likelihood weights.
/// @return Negative local loglikelihood (scalar). 
template<class Type>
Type GumbelNLL(vector<Type> lu1, vector<Type> lu2,
	       vector<Type> xc, vector<Type> beta, vector<Type> wgt) {
  // parameter on regular scale
  vector<Type> theta = 1.0 + (beta[0] + beta[1] * xc).exp();
  // Type pen = exp(beta[0]); // penalty term
  // pen = .005 * (exp(.01/pen) - 1.0); 
  // pre-computations
  vector<Type> logA = ((theta * lu1).exp() + (theta * lu2).exp()).log();
  vector<Type> thA = (logA/theta).exp();
  // loglikelihood
  vector<Type> ll = (theta-1.0) * (lu1 + lu2) - thA + (2.0/theta - 2.0) * logA;
  ll += (1 + (theta - 1.0)/thA).log();
  ll *= wgt;
  return -ll.sum(); // + pen;
}
