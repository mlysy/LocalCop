/// @file student.hpp
///
/// @brief A simple implementation of the Student-t CDF and quantile functions.
///
/// These are defined in terms of the incomplete beta function: <https://en.wikipedia.org/wiki/Student%27s_t-distribution#Cumulative_distribution_function>.
///
/// @todo Make these identical to their R implementations, which handle a few edge cases differently.

#ifndef LOCALCOP_STUDENT_HPP
#define LOCALCOP_STUDENT_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {

  /// Distribution function of the Student-t distribution.
  ///
  /// @param df Degrees of freedom.
  template <class Type>
  Type pt(Type q, Type df) {
    Type res = Type(0.5) * pbeta(df/(q*q + df), Type(0.5) * df, Type(0.5));
    return CppAD::CondExpLt(q, Type(0.0), res, Type(1.0) - res); 
  }
  VECTORIZE2_tt(pt)

  /// Quantile function of the Student-t distribution. 
  ///
  /// @param df Degrees of freedom.
  template <class Type>
  Type qt(Type p, Type df) {
    Type p2 = CppAD::CondExpGe(p, Type(0.5), p, Type(1.0) - p);
    Type res = qbeta(Type(2.0)  * (Type(1.0) - p2), Type(0.5) * df, Type(0.5));
    res = sqrt(df/res - df);
    return CppAD::CondExpGe(p, Type(0.5), res, -res);
  }
  VECTORIZE2_tt(qt)

  
} // end namespace LocalCop

#endif
