/// @file nll.hpp

#ifndef LOCALCOP_NLL_HPP
#define LOCALCOP_NLL_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop { 

/// Calculate negative log-likelihood calculations with censored data
///
/// Computes the negative of
///
/// ```
/// sum_{i in uncensored} w[i] * log_dCopula(u1[i], u2[i], X[i,] * eta) + 
/// sum_{i in first censored} w[i] * log_hCopula(u1[i], u2[i], X[i,] * eta) + 
/// sum_{i in second censored} ... + 
/// sum_{i in both censored} w[i] * log_pCopula(u1[i], u2[i], X[i,] * eta)
/// ```
///
/// @param[in] u1 First uniform variable.
/// @param[in] u2 Second uniform variable.
/// @param[in] w Weight vector.
/// @param[in] X Covariate matrix.
/// @param[in] eta Calibration parameter matrix.
/// @param[in] status_start Vector of length 4 giving the start locations for
/// uncensored, first censored, second censored, and both censored.
/// @param[in] status_length Same but giving length of each. 
/// @param Copula Object containing the d/p/h function for the given `Family` copula,
/// in the calibration scale (eta).
///
/// @return Value of the negative log-likelihood.
///
/// @todo Is it worth preallocating memory for the linear algebra operations?
template <class Type, class Family>
Type nll(cRefVector<Type>& u1, 
         cRefVector<Type>& u2, 
         cRefVector<Type>& w,
         cRefMatrix<Type>& X,
         cRefMatrix<Type>& eta,
         const std::vector<int>& status_start,
         const std::vector<int>& status_length
         Family<Type>& Copula) {
  Type res = Type(0.0);
  // uncensored
  int j=0;
  if(status_length(j) > 0) {
    res += (w.segment(status_start(j), status_length(j)) * Copula.dfun(
      u1.segment(status_start(j), status_length(j)),
      u2.segment(status_start(j), status_length(j)),
      X.block(status_start(j), 0, status_length(j), X.cols()) * eta,
      1 // give_log
    )).sum();
  }
  // first censored
  int j=1;
  if(status_length(j) > 0) {
    res += (w.segment(status_start(j), status_length(j)) * Copula.hfun(
      u1.segment(status_start(j), status_length(j)),
      u2.segment(status_start(j), status_length(j)),
      X.block(status_start(j), 0, status_length(j), X.cols()) * eta,
      1 // give_log
    )).sum();
  }
  // second censored
  int j=2;
  if(status_length(j) > 0) {
    res += (w.segment(status_start(j), status_length(j)) * Copula.hfun(
      u2.segment(status_start(j), status_length(j)), // switched u1 and u2
      u1.segment(status_start(j), status_length(j)),
      X.block(status_start(j), 0, status_length(j), X.cols()) * eta,
      1 // give_log
    )).sum();
  }
  // both censored
  int j=3;
  if(status_length(j) > 0) {
    res += (w.segment(status_start(j), status_length(j)) * Copula.pfun(
      u1.segment(status_start(j), status_length(j)),
      u2.segment(status_start(j), status_length(j)),
      X.block(status_start(j), 0, status_length(j), X.cols()) * eta,
      1 // give_log
    )).sum();
  }
  return -res;
}
}  
