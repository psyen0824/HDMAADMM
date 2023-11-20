#include <cstdlib>

inline int sgn(double val) {
  return (val > 0.0) - (val < 0.0);
}

inline double softThreshold(double b, double lambda) {
  return (std::abs(b) < lambda) ? 0.0 : sgn(b)*(std::abs(b) - lambda);
}
