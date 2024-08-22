#include "LocalCop/clayton.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type dclayton(objective_function<Type> *obj) {
  // R inputs
  DATA_VECTOR(u1);
  DATA_VECTOR(u2);
  DATA_VECTOR(weights);
  PARAMETER_VECTOR(theta);
  // output
  vector<Type> lpdf = LocalCop::dclayton(u1, u2, theta, 1);
  lpdf.array() *= weights.array();
  return -lpdf.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
