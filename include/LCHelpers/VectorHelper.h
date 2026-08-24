#ifndef LC_CONTENT_VECTOR_HELPER_H
#define LC_CONTENT_VECTOR_HELPER_H 1

#include "Objects/CartesianVector.h"
#include <cmath>

namespace lc_content {
namespace VectorHelper {
  inline float deltaPhi(float phi1, float phi2) {
    float dphi = phi1 - phi2;
    
    // wrap into [-pi, pi]
    if (dphi >  static_cast<float>(M_PI))
      dphi -= 2.f * static_cast<float>(M_PI);
    if (dphi < -static_cast<float>(M_PI))
      dphi += 2.f * static_cast<float>(M_PI);
    return dphi;
  } // deltaPhi

  inline float AngularDR2(float t1, float p1, float t2, float p2) {
    const float dt = t1 - t2;
    float dp = deltaPhi(p1, p2);

    return dt*dt + dp*dp;
  } // AngularDR2

  inline float AngularDR(float t1, float p1, float t2, float p2) {
    return std::sqrt(AngularDR2(t1, p1, t2, p2));
  } // AngularDR
} // namespace VectorHelper
} // namespace lc_content

#endif // #ifndef LC_CONTENT_VECTOR_HELPER_H