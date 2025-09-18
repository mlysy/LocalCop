#include "LocalCop/pbvn.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type pbvn(objective_function<Type> *obj)
{
  PARAMETER_VECTOR(b1);
  PARAMETER_VECTOR(b2);
  PARAMETER_VECTOR(rho);

  Type tiny = 0.0; // Set to 1e-12 for robustness
  vector<Type> ans = LocalCop::pbvn(b1, b2, rho) + tiny;
  return ans.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
