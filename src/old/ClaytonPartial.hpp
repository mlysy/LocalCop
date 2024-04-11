#include "LocalCop/clayton.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
  Type ClaytonPartial(objective_function<Type> *obj) {
    // R inputs
    DATA_VECTOR(u1);
    DATA_VECTOR(u2);
    DATA_VECTOR(weights)
    PARAMETER_VECTOR(theta);
    // output
    vector<Type> lpart(u1.size());
    LocalCop::Clayton<Type> cop;
    cop.set_u(u1, u2);
    cop.logpart(lpart, theta);
    lpart.array() *= weights.array();
    REPORT(lpart);
    return -lpart.sum();
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this