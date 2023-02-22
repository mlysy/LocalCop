

#include "LocalCop/clayton.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type ClaytonNLL(objective_function<Type> *obj) {
  // R inputs
  DATA_VECTOR(u1);
  DATA_VECTOR(u2);
  DATA_VECTOR(weights)
  PARAMETER_VECTOR(theta);
  // output
  vector<Type> lpdf(u1.size());
  LocalCop::Clayton<Type> cop;
  cop.set_u(u1, u2);
  cop.logpdf(lpdf, theta);
  lpdf.array() *= weights.array();
  return -lpdf.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
