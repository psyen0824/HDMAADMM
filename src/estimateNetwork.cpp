#include <RcppEigen.h>

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

inline double softThreshold(double b, double lambda) {
  return (abs(b) < lambda) ? 0.0 : sgn(b)*(abs(b) - lambda);
}

// [[Rcpp::export]]
Rcpp::List estimateNetwork(
    Eigen::Map<Eigen::MatrixXd> X,
    Eigen::Map<Eigen::MatrixXd> Y,
    Eigen::Map<Eigen::MatrixXd> M1,
    Eigen::Map<Eigen::MatrixXd> alpha,
    Eigen::Map<Eigen::MatrixXd> beta,
    Eigen::Map<Eigen::MatrixXd> gamma,
    Eigen::Map<Eigen::MatrixXd> tauAlpha,
    Eigen::Map<Eigen::MatrixXd> tauBeta,
    double rho,
    double lambda1a,
    double lambda1b,
    double lambda1g,
    double lambda2a,
    double lambda2b,
    Eigen::Map<Eigen::MatrixXd> XtX,
    Eigen::Map<Eigen::MatrixXd> XtM1,
    Eigen::Map<Eigen::MatrixXd> M1tM1PlusRhoInv,
    Eigen::Map<Eigen::MatrixXd> M1tY,
    Eigen::Map<Eigen::MatrixXd> XtY,
    Eigen::Map<Eigen::MatrixXd> L
) {
  int p = M1.cols(), p2 = X.cols(), j;
  Eigen::MatrixXd alphaStep1 = (XtX + Eigen::MatrixXd::Identity(p2, p2)).llt().solve(XtM1 + rho*alpha - tauAlpha);
  Eigen::MatrixXd betaStep2 = M1tM1PlusRhoInv * (M1tY - XtM1.transpose() * gamma + rho*beta - tauBeta);

  double numerator, crossProd;
  Eigen::MatrixXd alphaNew = alpha, betaNew = beta;
  for (j = 0; j < p; j++) {
    crossProd = alphaNew.row(0).dot(L.col(j)) - alphaNew(0, j) * L(j, j);
    numerator = softThreshold(lambda2a * crossProd + tauAlpha(0, j) + rho * alphaStep1(0, j), lambda1a);
    alphaNew(0, j) = numerator / (lambda2a * (1 - L(j, j)) + rho);
  }

  for (j = 0; j < p; j++) {
    crossProd = betaNew.col(0).dot(L.col(j)) - betaNew(j, 0) * L(j, j);
    numerator = softThreshold(lambda2b * crossProd + tauBeta(j, 0) + rho * betaStep2(j, 0), lambda1b);
    betaNew(j, 0) =  numerator / (lambda2b * (1 - L(j, j)) + rho);
  }

  Eigen::MatrixXd gammaTemp = XtY - XtM1 * betaStep2;
  for (j = 0; j < p2; j++) {
    gammaTemp(j, 0) = softThreshold(gammaTemp(j, 0), lambda1g);
  }
  Eigen::MatrixXd gammaNew = XtX.llt().solve(gammaTemp);
  Eigen::MatrixXd tauAlphaNew = tauAlpha + rho * (alphaStep1-alphaNew);
  Eigen::MatrixXd tauBetaNew = tauBeta + rho * (betaStep2-betaNew);

  return Rcpp::List::create(
    Rcpp::Named("alphaStep1") = alphaStep1,
    Rcpp::Named("betaStep2") = betaStep2,
    Rcpp::Named("alpha") = alphaNew,
    Rcpp::Named("beta") = betaNew,
    Rcpp::Named("gamma") = gammaNew,
    Rcpp::Named("tauAlpha") = tauAlphaNew,
    Rcpp::Named("tauBeta") = tauBetaNew
  );
}
