#include "LocalCop/student.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type qt(objective_function<Type> *obj) {
  // R inputs
  PARAMETER_VECTOR(p);
  PARAMETER_VECTOR(df);
  // output
  vector<Type> res = LocalCop::qt(p, df);
  ADREPORT(res);
  return sum(res);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
