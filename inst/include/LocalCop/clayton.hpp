/// @file clayton.hpp

#ifndef LOCALCOP_CLAYTON_HPP
#define LOCALCOP_CLAYTON_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {
  
  /// Calculate Clayton copula CDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Clayton copula with the range $[0,\infty)$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula CDF.
  template <class Type>
  Type pclayton(Type u1, Type u2, Type theta, int give_log=0) {
    Type logans = -Type(1.0)/theta * log(pow(u1, -theta) + pow(u2, -theta) - Type(1.0));
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(pclayton)    
    
    
  /// Calculate Clayton copula partial derivative with respect to u1.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Clayton copula with the range $[0,\infty)$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the h-function.
  template <class Type>
  Type hclayton(Type u1, Type u2, Type theta, int give_log=0) {
    Type logans = -(Type(1.0)+theta) * log(u1);
    logans -= (Type(1.0)+Type(1.0)/theta) * log(pow(u1, -theta) + pow(u2, -theta) - Type(1.0));
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(hclayton)    
      
  /// Calculate Clayton copula PDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Clayton copula with the range $[0,\infty)$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula PDF.
  template <class Type>
  Type dclayton(Type u1, Type u2, Type theta, int give_log=0) {
    Type logans = log(Type(1.0) + theta) - (Type(1.0) + theta) * (log(u1) + log(u2));
    logans -= (Type(2.0) + Type(1.0)/theta) * log(pow(u1, -theta) + pow(u2, -theta) - Type(1.0));
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(dclayton)
      
} // end namespace LocalCop

#endif // LOCALCOP_CLAYTON_HPP
