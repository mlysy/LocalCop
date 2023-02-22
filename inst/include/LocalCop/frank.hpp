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
    Type ans = log(Type(1.0) + (exp(-theta*u1) - Type(1.0)) * (exp(-theta*u2) - Type(1.0)) / (exp(-theta) - Type(1.0)));
    ans *= - Type(1.0)/theta;
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(pfrank)


  /// Frank Copula model.
  ///
  /// Usage:
  /// ```
  /// Frank<Type> cop; // create `cop` object
  /// cop.logpdf(lpdf, u1, u2, theta); // use `logpdf()` method of `cop`
  /// ```
  template <class Type>
  class Frank {
  private:
    // internal memory
    Vector_t<Type> u1_; // u1 and u2
    Vector_t<Type> u2_;
    Vector_t<Type> num_; // for internal calculations
    Vector_t<Type> den_;
    void exp_prod(RefVector_t<Type> ep, cRefVector_t<Type>& theta);
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

  /// Calculate `(1 - exp(-theta*u1)) * (1 - exp(-theta*u2))`.
  ///
  /// @param[out] ep Output vector.
  /// @param[in] theta Parameter vector.
  ///
  /// @warning Requires internal values of `u1_` and `u2_` to be set.
  template <class Type>
  inline void Frank<Type>::exp_prod(RefVector_t<Type> ep, 
                                    cRefVector_t<Type>& theta) {
    ep = (Type(1.0) - (-theta.array() * u1_.array()).exp()) * 
      (Type(1.0) - (-theta.array() * u2_.array()).exp());
    return;
  }
  
  /// Set the internal values of `u1` and `u2`.
  ///
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  template <class Type>
  inline void Frank<Type>::set_u(cRefVector_t<Type>& u1,
				   cRefVector_t<Type>& u2) {
    u1_ = u1.array();
    u2_ = u2.array();
    return;
  }


  /// @param[out] lpdf Vector of calculated log-densities.
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Frank<Type>::logpdf(RefVector_t<Type> lpdf, 
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
  inline void Frank<Type>::logpdf(RefVector_t<Type> lpdf, 
				    cRefVector_t<Type>& theta) {
    num_ = Type(1.0) - theta.array().exp();
    exp_prod(den_, theta);
    den_ = num_ - den_;
    den_.array() *= den_.array();
    num_.array() *= theta.array();
    lpdf = (num_.array() / den_.array()).log() - theta.array() * (u1_ + u2_).array();
    return;
  }

  /// @param[out] lcdf Vector of calculated log-CDFs.
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Frank<Type>::logcdf(RefVector_t<Type> lcdf, 
				    cRefVector_t<Type>& u1,
				    cRefVector_t<Type>& u2,
				    cRefVector_t<Type>& theta) {
    set_u(u1, u2);
    logcdf(lcdf, theta);
  }

  /// @param[out] lcdf Vector of calculated log-CDFs.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Frank<Type>::logcdf(RefVector_t<Type> lcdf, 
				    cRefVector_t<Type>& theta) {
    exp_prod(num_, theta);
    den_ =  theta.array().exp() - Type(1.0);
    lcdf = (Type(1.0) + (num_.array() / den_.array())).log();
    lcdf.array() *= (-Type(1.0)/theta.array());
    lcdf = lcdf.array().log();
    return;
  }
  

  
  /// @param[out] lpart Vector of calculated log-partial copula derivatives.
  /// @param[in] u1 First uniform response vector.
  /// @param[in] u2 Second uniform response vector.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Frank<Type>::logpart(RefVector_t<Type> lpart, 
                                    cRefVector_t<Type>& u1,
                                    cRefVector_t<Type>& u2,
                                    cRefVector_t<Type>& theta) {
    set_u(u1, u2);
    logpart(lpart, theta);
  }
  
  /// @param[out] lpart Vector of calculated log-partial copula derivatives.
  /// @param[in] theta Vector of copula parameters.
  template <class Type>
  inline void Frank<Type>::logpart(RefVector_t<Type> lpart, 
                                   cRefVector_t<Type>& theta) {
    exp_prod(den_, theta);
    den_ = theta.array().exp() -Type(1.0) + den_.array();
    num_ = (-theta.array() * u1_.array()).exp();
    num_.array() *= (-theta.array() * u2_.array()).exp() - Type(1.0);
    lpart = (num_.array() / den_.array()).log();
    return;
  }
  
  
  
  

} // end namespace LocalCop

#endif
