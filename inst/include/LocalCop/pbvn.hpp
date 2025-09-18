/// @file pbvn.hpp

#ifndef LOCALCOP_PBVN_HPP
#define LOCALCOP_PBVN_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

namespace LocalCop {

  /*
    This class evaluates the integral of the bivariate normal distribution via the direct method.
    Since we already have pnorm which is fast and accurate, we can avoid doing 2D integration by
    conditioning on one of the two variables, in our case on X1 and its associated variable b1.
  */
  template<class Float>
  struct BVNIntegrand {
    typedef Float Scalar; // Required by integrate
    Float b1, b2, rho;         // Parameters 
    // Evaluate conditional CDF
    Float operator() (Float x) {
      Float loc = rho * x;
      Float scale = sqrt(1 - rho * rho);

      // Replace Float(0.0) by adding mu to parameters above to control mean
      // Replace the first Float(1.0) by adding sigma1 to parameters above to control the first s.d.
      // Replace the second Float(1.0) by adding sigma2 to parameters above to control the second s.d.
      Float ans = pnorm((b2 - loc) / scale, Float(0.0), Float(1.0)) *
                  dnorm(x, Float(0.0), Float(1.0), false);


      // Float ans = pnorm((b2 - loc) / scale, Float(0.0), Float(1.0), false, true) *
      //             dnorm(x, Float(0.0), Float(1.0), false);
      return ans;
    }
    // Integrate conditional CDF into the CDF
    Float integrate() {
      using gauss_kronrod::integrate;
      Float ans =
        integrate(*this, Float(-INFINITY), b1);
      return ans;
    }
  };


  // An externally available integration evaluator
  template<class Float>
  Float pbvn(Float b1, Float b2, Float rho) {
    BVNIntegrand<Float> f = {b1, b2, rho};
    return f.integrate();
  }

  VECTORIZE3_ttt(pbvn)
}

#endif // LOCALCOP_PBVN_HPP
