#include "LocalCop/student.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type dstudent(objective_function<Type> *obj) {
  // R inputs
  DATA_VECTOR(u1);
  DATA_VECTOR(u2);
  DATA_VECTOR(weights);
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(nu);
  // output
  vector<Type> lpdf = LocalCop::dstudent(u1, u2, theta, nu, 1);
  lpdf.array() *= weights.array();
  return -lpdf.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
