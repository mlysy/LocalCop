/// @file gumbel.hpp

#ifndef LOCALCOP_GUMBEL_HPP
#define LOCALCOP_GUMBEL_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {

  /// Calculate Gumbel copula CDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Gumbel copula with the range $[1,\infty)$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula CDF.
  template <class Type>
  Type pgumbel(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 =  pow(-log(u1), theta) ;
    Type term2 =  pow(-log(u2), theta) ;
    Type logans = -pow( term1 + term2, Type(1.0)/theta);
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(pgumbel)

  /// Calculate Gumbel copula partial derivative wrt u1.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Gumbel copula with the range $[1,\infty)$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the h-function.
  template <class Type>
  Type hgumbel(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 = pow(-log(u1), theta) + pow(-log(u2), theta);
    Type term2 = -log(u1);
    Type term3 = Type(1.0) / theta;
    Type logans = -pow(term1, term3) + (theta - Type(1.0)) * log(term2) + (term3 - Type(1.0)) * log(term1) + term2;
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(hgumbel)
      
  /// Calculate Gumbel copula PDF.
  //
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Gumbel copula with the range $[1,\infty)$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula PDF.
  template <class Type>
  Type dgumbel(Type u1, Type u2, Type theta, int give_log=0) {
    Type log_u1 = log(u1);
    Type log_u2 = log(u2);
    Type loglog_u1 = log(-log_u1);
    Type loglog_u2 = log(-log_u2);
    Type log_theta = log(theta - 1.0);
    Type lsum = logspace_add(theta * loglog_u1, theta * loglog_u2);
    Type ans = (theta - 1.0) * (loglog_u1 + loglog_u2);
    ans += (2.0 * (1.0/theta - 1.0)) * lsum - exp(1.0/theta * lsum);
    ans += log_theta + logspace_add(-log_theta,  -1.0/theta * lsum);
    ans -= log_u1 + log_u2;
    if(give_log) return ans; else return exp(ans);
  }
  VECTORIZE4_ttti(dgumbel)

} // end namespace LocalCop

#endif // LOCALCOP_GUMBEL_HPP
