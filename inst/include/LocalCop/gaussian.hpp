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
  template<class Float>
  struct gaussianIntegrate {
    typedef Float Scalar;
    Float theta;
    Float u1, u2;

    Float operator() () {
      Float ans = dgaussian(u1, u2, theta, 0);
      return ans;
    }

    Float integrate(Float x1, Float x2) {
      using gauss_kronrod::mvIntegrate;
      Float ans = mvIntegrate(*this).
        wrt(u1, -INFINITY, x1).
        wrt(u2, -INFINITY, x2) ();
      return ans;
    }
  };

  template <class Type>
  Type pgaussian_eval(Type u1, Type u2, Type theta) {
    gaussianIntegrate<Type> gi;
    gi.theta = theta;
    Type ans = gi.integrate(qnorm(u1), qnorm(u2));
    return ans;
  }

  TMB_BIND_ATOMIC(pgaussian_func, 001, pgaussian_eval(x[0], x[1], x[2]))

  /// Calculate Gaussian copula CDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Gaussian copula with the range $(-1, 1)$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the p-function.  
  template<class Type>
  Type pgaussian(Type u1, Type u2, Type theta, int give_log=0) {
    vector<Type> args(4); // Last index reserved for derivative order
    args << u1, u2, theta, 0;
    Type ans = Type(LocalCop::pgaussian_func(CppAD::vector<Type>(args))[0]); // TODO: edited
    if(give_log) return log(ans); else return ans;
  }

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
    Type z1 = qnorm(u1);
    Type z2 = qnorm(u2);
    Type determinant = Type(1.0) - theta * theta;
    Type ans = pnorm((z2 - theta * z1) / sqrt(determinant));
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(hgaussian)
      
  /// Calculate Gaussian copula PDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Gaussian copula with the range $(-1, 1)$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the copula PDF. 
  template <class Type>
  Type dgaussian(Type u1, Type u2, Type theta, int give_log=0) {
    // normal quantiles
    Type z1 = qnorm(u1);
    Type z2 = qnorm(u2);
    Type det = 1.0 - theta*theta;
    Type ans = theta*theta * (z1*z1 + z2*z2) - 2.0*theta * z1*z2;
    ans = -.5 * (ans / det + log(det));
    if(give_log) return ans; else return exp(ans);
  }
  VECTORIZE4_ttti(dgaussian)

} // end namespace LocalCop

#endif // LOCALCOP_GAUSSIAN_HPP
