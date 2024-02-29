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
  /// @param[in] theta Parameter of the Gumbel copula with the range $[1,\infty]$.
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
  /// @param[in] theta Parameter of the Gumbel copula with the range $[1,\infty]$.
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
  /// @param[in] theta Parameter of the Gumbel copula with the range $[1,\infty]$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula PDF.
  template <class Type>
  Type dgumbel(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 = log(u1);
    Type term2 = log(u2);
    Type term3 = pow(-term1, theta) + pow(-term2, theta);
    Type term4 = Type(1.0) / theta;
    Type term5 = theta - Type(1.0);
    Type term6 = term4 - Type(1.0);
    Type logans = term5 * (log(term1) + log(term1)) + Type(2.0) * term6 * log(term3) - pow(term3, term4) - term1 - term2;
    logans += log(Type(1.0) + (theta - Type(1.0)) * pow(term3, -term4))
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(dgumbel)

} // end namespace LocalCop

#endif // LOCALCOP_GUMBEL_HPP
