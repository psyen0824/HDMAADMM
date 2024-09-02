#include <RcppEigen.h>
#include "utility.h"

// [[Rcpp::export]]
Rcpp::List getLogLikelihood(
    Eigen::MatrixXd X,
    Eigen::MatrixXd Y,
    Eigen::MatrixXd M1,
    Eigen::MatrixXd alpha,
    Eigen::MatrixXd beta,
    Eigen::MatrixXd gamma
) {
  Eigen::MatrixXd A = M1 - X * alpha;
  double l1 = (-0.5) * (A.transpose() * A).diagonal().array().sum();
  Eigen::MatrixXd B = Y - X * gamma - M1 * beta;
  Eigen::MatrixXd l2 = (-0.5) * B.transpose() * B;
  return Rcpp::List::create(
    Rcpp::Named("l") = l1+l2(0, 0),
    Rcpp::Named("l1") = l1,
    Rcpp::Named("l2") = l2(0, 0)
  );
}

double getObjective(
    Eigen::MatrixXd X,
    Eigen::MatrixXd Y,
    Eigen::MatrixXd M1,
    Eigen::MatrixXd alpha,
    Eigen::MatrixXd beta,
    Eigen::MatrixXd gamma,
    int penaltyType,
    double lambda1a,
    double lambda1b,
    double lambda1g,
    double lambda2a,
    double lambda2b,
    double kappa,
    Eigen::MatrixXd laplacianMatrixA,
    Eigen::MatrixXd laplacianMatrixB,
    double lambda2aStar,
    double lambda2bStar
) {
  Rcpp::List logLikelihoodList = getLogLikelihood(X, Y, M1, alpha, beta, gamma);
  double logLik = logLikelihoodList("l");

  double p1 = lambda1g * gamma.array().abs().sum();
  double p2 = lambda1a * alpha.array().abs().sum() + lambda1b * beta.array().abs().sum();
  double p3 = 0.0;
  if (penaltyType == 1) {
    // Elastic Net
    p3 = lambda2a * alpha.row(0).dot(alpha.row(0)) + lambda2b * beta.col(0).dot(beta.col(0));
  } else if (penaltyType == 2) {
    // Network
    Eigen::MatrixXd alphaTemp = alpha * laplacianMatrixA * alpha.transpose();
    Eigen::MatrixXd betaTemp = beta.transpose() * laplacianMatrixB * beta;
    p3 = lambda2a * alphaTemp(0, 0) + lambda2b * betaTemp(0, 0);
  } else if (penaltyType == 3) {
    // Pathway Lasso
    p3 = (alpha.row(0).transpose().array() * beta.col(0).array()).abs().sum() +
      lambda2a * alpha.row(0).dot(alpha.row(0)) + lambda2b * beta.col(0).dot(beta.col(0));
    p3 *= kappa;
  }  else if (penaltyType == 4) {
    // Pathway Network
    Eigen::MatrixXd alphaTemp = alpha * laplacianMatrixA * alpha.transpose();
    Eigen::MatrixXd betaTemp = beta.transpose() * laplacianMatrixB * beta;
    p3 = (alpha.row(0).transpose().array() * beta.col(0).array()).abs().sum() +
      lambda2a * alphaTemp(0, 0) + lambda2a * lambda2aStar * alpha.row(0).dot(alpha.row(0)) +
      lambda2b * betaTemp(0, 0) + lambda2b * lambda2bStar * beta.col(0).dot(beta.col(0));
    p3 *= kappa;
  }
  return logLik + p1 + p2 + p3;
}
