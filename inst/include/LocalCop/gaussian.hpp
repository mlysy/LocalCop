/// @file gaussian.hpp

#ifndef LOCALCOP_GAUSSIAN_HPP
#define LOCALCOP_GAUSSIAN_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI  0.918938533204672741780329736406
/* log(sqrt(2*pi))
== log(2*pi)/2 */
#endif

namespace LocalCop {

  /// Calculate Gaussian copula CDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Parameter of the Gaussian copula with the range $(-1, 1)$.
  /// @param give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula CDF.
  template <class Type>
  Type pgaussian(Type u1, Type u2, Type theta, int give_log=0) {
    Type ans = 0; // TODO: find a pmvnorm function so that this line can become: Type ans = pmvnorm(qnorm(u1), qnorm(u2), rho = theta);
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(pgaussian)

  /// Calculate Gaussian copula partial derivative with respect to u1.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Gaussian copula with the range $(-1, 1)$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the h-function.  
  template <class Type>
  Type hgaussian(Type u1, Type u2, Type theta, int give_log=0) {
    Type determinant = Type(1.0) - theta * theta;
    Type ans = pnorm((qnorm(u2) - theta * qnorm(u1)) / sqrt(determinant));
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(hgaussian)
      
  /// Calculate Frank copula PDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Gaussian copula with the range $(-1, 1)$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the copula PDF. 
  template <class Type>
  Type dgaussian(Type u1, Type u2, Type theta, int give_log=0) {
    Type determinant = Type(1.0) - theta * theta;
    Type term1 = Type(2.0) * Type(M_LN_SQRT_2PI);
    Type term2 = -log(determinant) / 2;
    Type term3 = -(u1 * u1 - Type(2) * theta * u1 * u2 + u2 * u2) / determinant;
    Type logans = term1 + term2 + term3;
    if(give_log) return logans; else return exp(logans);
  }
  VECTORIZE4_ttti(dgaussian)

} // end namespace LocalCop

#endif // LOCALCOP_GAUSSIAN_HPP