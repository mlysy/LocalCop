/// @file frank.hpp

#ifndef LOCALCOP_FRANK_HPP
#define LOCALCOP_FRANK_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {

  /// Calculate Frank copula CDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Frank copula with the range $R \setminus \{0\}$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula CDF.
  template <class Type>
  Type pfrank(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 =  exp(-theta*u1) - Type(1.0) ;
    Type term2 =  exp(-theta*u2) - Type(1.0) ; 
    Type term3 =  exp(-theta) - Type(1.0) ;
    Type ans = log(Type(1.0) + (term1 * term2) / term3);
    ans /= Type(-1.0) * theta;
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(pfrank)

  /// Calculate Frank copula partial derivative with respect to u1.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Frank copula with the range $R \setminus \{0\}$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the h-function.  
  template <class Type>
  Type hfrank(Type u1, Type u2, Type theta, int give_log=0) {
    Type exp_u2 = exp(-theta*u2) - 1.0;
    Type ans = exp(-theta*u1) * exp_u2;
    ans /= ((exp(-theta) - 1.0) + ans - exp_u2);
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(hfrank)
      
  /// Calculate Frank copula PDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Frank copula with the range $R \setminus \{0\}$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the copula PDF. 
  template <class Type>
  Type dfrank(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 =  exp(-theta*u1);
    Type term2 =  exp(-theta*u2); 
    Type term3 =  exp(-theta) - 1.0;
    Type term4 = (term3 + (term1 - 1.0) * (term2 - 1.0));
    term4 *= term4;
    Type ans = -theta * term3 * term1 * term2;
    ans /= term4;
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(dfrank)

} // end namespace LocalCop

#endif // LOCALCOP_FRANK_HPP
