/// @file clayton.hpp

#ifndef LOCALCOP_CLAYTON_HPP
#define LOCALCOP_CLAYTON_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {

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
      
      
  /// Clayton Copula model.
  ///
  /// Usage:
  /// ```
  /// Clayton<Type> cop; // create `cop` object
  /// cop.logpdf(lpdf, u1, u2, theta); // use `logpdf()` method of `cop`
  /// ```
  template <class Type>
  class Clayton {
  private:
    // internal memory
    Vector_t<Type> lu1_; // log of u1 and u2
    Vector_t<Type> lu2_;
  public:
    // Set response vectors.
    void set_u(cRefVector_t<Type>& u1,
               cRefVector_t<Type>& u2);
    // Calculate log-density for vector inputs.
    void logpdf(RefVector_t<Type> lpdf, 
                cRefVector_t<Type>& u1,
                cRefVector_t<Type>& u2,
                cRefVector_t<Type>& theta);
    // Calculate log-density for pre-allocated response vectors.
    void logpdf(RefVector_t<Type> lpdf, 
                cRefVector_t<Type>& theta);
    // Calculate log-CDF for vector inputs.
    void logcdf(RefVector_t<Type> lcdf, 
                cRefVector_t<Type>& u1,
                cRefVector_t<Type>& u2,
                cRefVector_t<Type>& theta);
    // Calculate log-CDF for pre-allocated response vectors.
    void logcdf(RefVector_t<Type> lcdf, 
                cRefVector_t<Type>& theta);
    // Calculate log-Partial for vector inputs.
    void logpart(RefVector_t<Type> lpart, 
                 cRefVector_t<Type>& u1,
                 cRefVector_t<Type>& u2,
                 cRefVector_t<Type>& theta);
    // Calculate log-Partial for pre-allocated response vectors.
    void logpart(RefVector_t<Type> lpart, 
                 cRefVector_t<Type>& theta);
  };
      
  /// Set the internal values of `lu1 = log(u1)` and `lu2 = log(u2)`.
  ///
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  template <class Type>
  inline void Clayton<Type>::set_u(cRefVector_t<Type>& u1,
                                   cRefVector_t<Type>& u2) {
    lu1_ = u1.array().log();
    lu2_ = u2.array().log();
    return;
  }
      
      
  /// @param[out] lpdf Vector of calculated log-densities.    
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Clayton<Type>::logpdf(RefVector_t<Type> lpdf, 
                                    cRefVector_t<Type>& u1,
                                    cRefVector_t<Type>& u2,
                                    cRefVector_t<Type>& theta) {
    set_u(u1, u2);
    logpdf(lpdf, theta);
    return;
  }
      
  /// @param[out] lpdf Vector of calculated log-densities.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Clayton<Type>::logpdf(RefVector_t<Type> lpdf, 
                                    cRefVector_t<Type>& theta) {
    lpdf = (Type(1.0) + theta.array()).log() - (Type(1.0) + theta.array()) * (lu1_ + lu2_).array();
    lpdf.array() -= (Type(2.0) + Type(1.0)/theta.array()) * ((-theta.array() * lu1_.array()).exp() + (-theta.array() * lu2_.array()).exp() - Type(1.0)).log();
    return;
  }
      
  /// @param[out] lcdf Vector of calculated log-CDFs.
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Clayton<Type>::logcdf(RefVector_t<Type> lcdf, 
                                    cRefVector_t<Type>& u1,
                                    cRefVector_t<Type>& u2,
                                    cRefVector_t<Type>& theta) {
    set_u(u1, u2);
    logcdf(lcdf, theta);
  }
      
  /// @param[out] lcdf Vector of calculated log-CDFs.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Clayton<Type>::logcdf(RefVector_t<Type> lcdf, 
                                    cRefVector_t<Type>& theta) {
    lcdf = (- theta.array() * lu1_.array()).exp() + (- theta.array() * lu2_.array()).exp() - Type(1.0);
    lcdf.array() = -Type(1.0)/theta.array() * lcdf.array().log();
    return;      
  }
      
      
      
  /// @param[out] lpart Vector of calculated log-partial copula derivatives.
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Clayton<Type>::logpart(RefVector_t<Type> lpart, 
                                     cRefVector_t<Type>& u1,
                                     cRefVector_t<Type>& u2,
                                     cRefVector_t<Type>& theta) {
    set_u(u1, u2);
    logpart(lpart, theta);
  }
      
  /// @param[out] lpart Vector of calculated log-partial copula derivatives.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Clayton<Type>::logpart(RefVector_t<Type> lpart, 
                                     cRefVector_t<Type>& theta) {
    lpart = (- theta.array() - Type(1.0))* lu1_.array(); 
    lpart.array() -= (Type(1.0) + Type(1.0)/theta.array()) * ((-theta.array() * lu1_.array()).exp() + (-theta.array() * lu2_.array()).exp() - Type(1.0)).log();
    return;
  }
      
      
      
      
} // end namespace LocalCop

#endif
