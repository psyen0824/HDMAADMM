#include <Rcpp.h>

inline template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

inline double softThreshold(double b, double lambda) {
  return (abs(b) < lambda) ? 0.0 : sgn(b)*(abs(b) - lambda);
}
