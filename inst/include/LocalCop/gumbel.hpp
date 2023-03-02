/// @file gumbel.hpp

#ifndef LOCALCOP_GUMBEL_HPP
#define LOCALCOP_GUMBEL_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {

  /// Calculate Gumbel copula CDF.
  ///
  /// Usage:
  /// ```
  /// lcdf = pgumbel<Type>(u1, u2, theta, 1); // 
  /// ```

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
  /// Usage:
  /// ```
  /// lpart = hgumbel<Type>(u1, u2, theta, 1); // 
  /// ```
    
  template <class Type>
  Type hgumbel(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 =  pow(-log(u1), theta) ;
    Type term2 =  pow(-log(u2), theta) ;
    Type term3 =  exp(-pow( term1 + term2, Type(1.0)/theta)) ;
    Type ans = pow(term1 + term2, Type(1.0)/theta - Type(1.0)) ;
    ans *= pow(-log(u1), theta - Type(1.0)) * term3 ;
    ans  /= u1;
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(hgumbel)
      
  /// Calculate Gumbel copula PDF.
  ///
  /// Usage:
  /// ```
  /// lpdf = dgumbel<Type>(u1, u2, theta, 1); // 
  /// ```
    
  template <class Type>
  Type dgumbel(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 =  pow(-log(u1), theta) ;
    Type term2 =  pow(-log(u2), theta) ;
    Type term3 =  exp(-pow( term1 + term2, Type(1.0)/theta)) ;
    Type ans = pow(log(u1)*log(u2), theta - Type(1.0)) ;   
    ans *= pow( term1 + term2, Type(2.0) * (Type(1.0)/theta - Type(1.0))) * term3 ;
    ans *= (Type(1.0) + (theta - Type(1.0)) * pow(term1 + term2, -Type(1.0)/theta)) ;
    ans /= (u1*u2) ;
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(dgumbel)    
    
  
    
  
  
  

} // end namespace LocalCop

#endif