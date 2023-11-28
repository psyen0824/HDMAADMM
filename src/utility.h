#include <RcppEigen.h>
#include <iostream>

inline int sgn(double val) {
  return (val > 0.0) - (val < 0.0);
}

inline double softThreshold(double b, double lambda) {
  return (std::abs(b) < lambda) ? 0.0 : sgn(b)*(std::abs(b) - lambda);
}

inline void printCoefficient(double *coef, std::string coefName, int numToPrint) {
  int i;
  for (i = 0; i < numToPrint; i++){
    Rcpp::Rcout << ", " <<  coefName << "[" << i+1 << "]: " << coef[i];
  }
}
