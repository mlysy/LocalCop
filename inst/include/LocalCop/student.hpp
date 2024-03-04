/// @file student.hpp

#ifndef LOCALCOP_STUDENT_HPP
#define LOCALCOP_STUDENT_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

namespace LocalCop {

  /// Distribution function of the Student-t distribution.
  ///
  /// This implementation is defined in terms of the incomplete beta function: <https://en.wikipedia.org/wiki/Student%27s_t-distribution#Cumulative_distribution_function>.
  ///
  /// @param[in] q Quantile.
  /// @param[in] df Degrees of freedom.
  ///
  /// @return Value of the CDF at `q`.
  ///
  /// @todo Make this function identical to its R implementation, which handles a few edge cases differently.
  template <class Type>
  Type pt(Type q, Type df) {
    Type res = Type(0.5) * pbeta(df/(q*q + df), Type(0.5) * df, Type(0.5));
    return CppAD::CondExpLt(q, Type(0.0), res, Type(1.0) - res); 
  }
  VECTORIZE2_tt(pt)

  /// Quantile function of the Student-t distribution. 
  ///
  /// Defined via direct inversion of `pt()`.
  /// 
  /// @param[in] p Probability.
  /// @param[in] df Degrees of freedom.
  ///
  /// @return Quantile of the Student-t corresponding to `p`.
  ///
  /// @todo Make this function identical to its R implementation.
  template <class Type>
  Type qt(Type p, Type df) {
    Type p2 = CppAD::CondExpGe(p, Type(0.5), p, Type(1.0) - p);
    Type res = qbeta(Type(2.0)  * (Type(1.0) - p2), Type(0.5) * df, Type(0.5));
    res = sqrt(df/res - df);
    return CppAD::CondExpGe(p, Type(0.5), res, -res);
  }
  VECTORIZE2_tt(qt)

  /// Calculate Student-t copula PDF.
  //
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable.
  /// @param[in] theta Correlation parameter of the Student-t copula with the range $(-1, 1)$.
  /// @param[in] nu Degrees of freedom parameter.
  /// @param[in] give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the copula PDF.
  template <class Type>
  Type dstudent(Type u1, Type u2, Type theta, Type nu, int give_log=0) {
    // student-t quantiles
    Type y1 = qt(u1, nu);
    Type y2 = qt(u2, nu);
    // bivariate student log-pdf
    Type det = 1.0 - theta*theta;
    Type ans = (y1*y1 + y2*y2 - 2.0*theta * y1*y2) / det;
    ans = -(.5 * nu + 1.0) * log(1.0 + ans/nu);
    ans -= 2.0 * M_LN_SQRT_2PI + .5 * log(det);
    // marginals
    ans -= dt(y1, nu, 1) + dt(y2, nu, 1);
    if(give_log) return ans; else return exp(ans);
  }
  VECTORIZE5_tttti(dstudent)

  /// Calculate Student-t copula partial derivative with respect to u1.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Correlation parameter of the Student-t copula with the range $(-1, 1)$.
  /// @param[in] nu Degrees of freedom parameter.
  /// @param[in] give_log Whether or not to return on the log scale.
  ///
  /// @return Value of the h-function.  
  template <class Type>
  Type hstudent(Type u1, Type u2, Type theta, Type nu, int give_log=0) {
    // student-t quantiles
    Type y1 = qt(u1, nu);
    Type y2 = qt(u2, nu);
    Type loc = theta * y2;
    Type det = 1.0 - theta*theta;
    Type nu1 = nu + 1.0;
    Type scale = sqrt((nu + y2*y2)/nu1 * det);
    Type z = (y1 - loc)/scale;
    Type ans = pt(z, nu1);
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE5_tttti(hstudent)

  
} // end namespace LocalCop

#endif
