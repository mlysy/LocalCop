/// @file frank.hpp

#ifndef LOCALCOP_FRANK_HPP
#define LOCALCOP_FRANK_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {

  /// Calculate Frank copula CDF.
  ///
  /// Usage:
  /// ```
  /// lcdf = pfrank<Type>(u1, u2, theta, 1); // 
  /// ```

  template <class Type>
  Type pfrank(Type u1, Type u2, Type theta, int give_log=0) {
    Type ans = log(Type(1.0) + ((exp(-theta*u1) - Type(1.0)) * (exp(-theta*u2)- Type(1.0))) / (exp(-theta)- Type(1.0)));
    ans /= Type(-1.0) * theta;
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(pfrank)

  /// Calculate Frank copula partial derivative wrt u1.
  ///
  /// Usage:
  /// ```
  /// lpart = hfrank<Type>(u1, u2, theta, 1); // 
  /// ```
    
  template <class Type>
  Type hfrank(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 =  exp(-theta*u1) - Type(1.0) ;
    Type term2 =  exp(-theta*u2) - Type(1.0) ; 
    Type term3 =  exp(-theta) - Type(1.0) ;
    Type ans = (term1 + Type(1.0)) * term2 ;
    ans  /= (term3 + term1 * term2);
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(hfrank)
      
  /// Calculate Frank copula PDF.
  ///
  /// Usage:
  /// ```
  /// lpdf = dfrank<Type>(u1, u2, theta, 1); // 
  /// ```
    
  template <class Type>
  Type dfrank(Type u1, Type u2, Type theta, int give_log=0) {
    Type term1 =  exp(-theta*u1) - Type(1.0) ;
    Type term2 =  exp(-theta*u2) - Type(1.0); 
    Type term3 =  exp(-theta) - Type(1.0) ;
    Type term4 = (term3 + term1 * term2) ;
    term4 *= term4 ;
    //std::cout << "term1 = " << term1 << ", term2 = " << term2 << std::endl;
    //std::cout << "term3 = " << term3 << ", term4 = " << term4 << std::endl;
    //printf("term1 = %f, term2 = %f, term3 = %f, term4 = %f\n", term1, term2, term3, term4);
    Type ans = -theta * term3 * (term1 + Type(1.0)) * (term2 + Type(1.0)) ;
    ans /= term4 ;
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(dfrank)    
    
  
    
  
  
  

} // end namespace LocalCop

#endif
