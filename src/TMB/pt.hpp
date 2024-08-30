#include "LocalCop/student.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type pt(objective_function<Type> *obj) {
  // R inputs
  PARAMETER_VECTOR(q);
  PARAMETER_VECTOR(df);
  // output
  vector<Type> res = LocalCop::pt(q, df);
  ADREPORT(res);
  return sum(res);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
