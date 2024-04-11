#ifndef LOCALCOP_CONFIG_HPP
#define LOCALCOP_CONFIG_HPP

#include <Eigen/Dense>

namespace LocalCop {
  /// Import namespaces inside namespace
  using namespace Eigen;

  /// Ref<Type> templates for Eigen function arguments.
  ///
  /// Basic idea:
  ///
  /// ```
  /// template <class Type>
  /// void foo(RefVector_t<Type> y, cRefMatrix_t<Type>& x);
  /// ```
  ///
  /// Here:
  ///
  /// - `y` is an input/output argument (so can be modified by the code) and is a vector (i.e., one column matrix).
  /// - `x` is an input argument (constant, so can't be modified by the code) and is a matrix.
  ///
  template <class Type>
  using Matrix_t = Matrix<Type, Dynamic, Dynamic>;
  template <class Type>
  using Vector_t = Matrix<Type, Dynamic, 1>;
  template <class Type>
  using RowVector_t = Matrix<Type, 1, Dynamic>;
  template <class Type>
  using RefMatrix_t = Ref <Matrix<Type, Dynamic, Dynamic> >;
  template <class Type>
  using cRefMatrix_t = const Ref <const Matrix<Type, Dynamic, Dynamic> >;
  template <class Type>
  using RefVector_t = Ref <Matrix<Type, Dynamic, 1> >;
  template <class Type>
  using cRefVector_t = const Ref <const Matrix<Type, Dynamic, 1> >;
  template <class Type>
  using RefRowVector_t = Ref <Matrix<Type, 1, Dynamic> >;
  template <class Type>
  using cRefRowVector_t = const Ref <const Matrix<Type, 1, Dynamic> >;

}

#endif
