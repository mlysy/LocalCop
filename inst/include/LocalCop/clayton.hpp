/// @file clayton.hpp

#ifndef LOCALCOP_CLAYTON_HPP
#define LOCALCOP_CLAYTON_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {
  
  /// Calculate Clayton copula CDF.
  ///
  /// Usage:
  /// ```
  /// lcdf = pclayton<Type>(u1, u2, theta, 1); // 
  /// ```
  template <class Type>
  Type pclayton(Type u1, Type u2, Type theta, int give_log=0) {
    Type logans = -Type(1.0)/theta * log(pow(u1, -theta) + pow(u2, -theta) - Type(1.0));
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(pclayton)    
    
    
  /// Calculate Clayton copula partial derivative wrt u1.
  ///
  /// Usage:
  /// ```
  /// lpart = hclayton<Type>(u1, u2, theta, 1); // 
  /// ```
  template <class Type>
  Type hclayton(Type u1, Type u2, Type theta, int give_log=0) {
    Type logans = -(Type(1.0)+theta) * log(u1);
    logans -= (Type(1.0)+Type(1.0)/theta) * log(pow(u1, -theta) + pow(u2, -theta) - Type(1.0));
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(hclayton)    
      
  /// Calculate Clayton copula PDF.
  ///
  /// Usage:
  /// ```
  /// lpdf = dclayton<Type>(u1, u2, theta, 1); // 
  /// ```
      
  template <class Type>
  Type dclayton(Type u1, Type u2, Type theta, int give_log=0) {
    Type logans = log(Type(1.0) + theta) - (Type(1.0) + theta) * (log(u1) + log(u2));
    logans -= (Type(2.0) + Type(1.0)/theta) * log(pow(u1, -theta) + pow(u2, -theta) - Type(1.0));
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(dclayton)
        
      
      
} // end namespace LocalCop

#endif
