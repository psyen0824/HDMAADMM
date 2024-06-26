#include <RcppEigen.h>
#include <tuple>
#include "utility.h"

Eigen::MatrixXd upadteAlphaElasticNet(
    Eigen::MatrixXd alphaStep1,
    Eigen::MatrixXd tauAlpha,
    double rho,
    double lambda1a,
    double lambda2a
) {
  int p = alphaStep1.cols(), j;
  Eigen::MatrixXd alphaNew(1, p);
  for (j = 0; j < p; ++j) {
    alphaNew(0, j) = softThreshold(tauAlpha(0, j) + rho*alphaStep1(0, j), lambda1a) / (lambda2a + rho);
  }
  return alphaNew;
};

Eigen::MatrixXd upadteBetaElasticNet(
    Eigen::MatrixXd betaStep2,
    Eigen::MatrixXd tauBeta,
    double rho,
    double lambda1b,
    double lambda2b
) {
  int p = betaStep2.rows(), j;
  Eigen::MatrixXd betaNew(p, 1);
  for (j = 0; j < p; ++j) {
    betaNew(j, 0) = softThreshold(tauBeta(j, 0) + rho*betaStep2(j, 0), lambda1b) / (lambda2b + rho);
  }
  return betaNew;
};

Eigen::MatrixXd upadteAlphaNetwork(
    Eigen::MatrixXd laplacianMatrixA,
    Eigen::MatrixXd alpha,
    Eigen::MatrixXd alphaStep1,
    Eigen::MatrixXd tauAlpha,
    double rho,
    double lambda1a,
    double lambda2a
) {
  int p = alpha.cols(), j;
  double numerator, crossProd;
  Eigen::MatrixXd alphaNew = alpha;
  for (j = 0; j < p; j++) {
    crossProd = alphaNew(0, j) * laplacianMatrixA(j, j) - alphaNew.row(0).dot(laplacianMatrixA.col(j));
    numerator = softThreshold(lambda2a * crossProd + tauAlpha(0, j) + rho * alphaStep1(0, j), lambda1a);
    alphaNew(0, j) = numerator / (lambda2a * laplacianMatrixA(j, j) + rho);
  }
  return alphaNew;
};

Eigen::MatrixXd upadteBetaNetwork(
    Eigen::MatrixXd laplacianMatrixB,
    Eigen::MatrixXd beta,
    Eigen::MatrixXd betaStep2,
    Eigen::MatrixXd tauBeta,
    double rho,
    double lambda1b,
    double lambda2b
) {
  int p = beta.rows(), j;
  double numerator, crossProd;
  Eigen::MatrixXd betaNew = beta;
  for (j = 0; j < p; ++j) {
    crossProd = betaNew(j, 0) * laplacianMatrixB(j, j) - betaNew.col(0).dot(laplacianMatrixB.col(j));
    numerator = softThreshold(lambda2b * crossProd + tauBeta(j, 0) + rho * betaStep2(j, 0), lambda1b);
    betaNew(j, 0) =  numerator / (lambda2b * laplacianMatrixB(j, j) + rho);
  }
  return betaNew;
};

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> upadteAlphaBetaPathwayLasso(
    Eigen::MatrixXd alphaStep1,
    Eigen::MatrixXd betaStep2,
    Eigen::MatrixXd tauAlpha,
    Eigen::MatrixXd tauBeta,
    double rho,
    double lambda1a,
    double lambda1b,
    double kappa,
    double lambda2a,
    double lambda2b
) {
  int p = alphaStep1.cols(), j;
  // To do a fair comparison, we introduce lambda2a, lambda2b to represent L2 norm penalty
  double phi1 = 2*kappa*lambda2a+rho, phi2 = 2*kappa*lambda2b+rho;
  // Below is the old code: using nu
  // double phi1 = 2*kappa*nu+rho, phi2 = 2*kappa*nu+rho;
  double muAlpha, muBeta, denominator, numeratorAlpha, numeratorBeta;
  Eigen::MatrixXd alphaNew(1, p), betaNew(p, 1);
  for (j = 0; j < p; ++j) {
    muAlpha = tauAlpha(0, j) + rho*alphaStep1(0, j);
    muBeta = tauBeta(j, 0) + rho*betaStep2(j, 0);

    if (kappa == 0.0) {
      alphaNew(0, j) = softThreshold(muAlpha, lambda1a) / phi1;
      betaNew(j, 0) = softThreshold(muBeta, lambda1b) / phi2;
    } else {
      denominator = phi1*phi2-kappa*kappa;

      numeratorAlpha = phi2*(muAlpha-lambda1a)-kappa*(muBeta-lambda1b);
      numeratorBeta = phi1*(muBeta-lambda1b)-kappa*(muAlpha-lambda1a);
      if ((numeratorAlpha > 0) && (numeratorBeta > 0)) {
        alphaNew(0, j) = numeratorAlpha / denominator;
        betaNew(j, 0) = numeratorBeta / denominator;
      } else {
        numeratorAlpha = phi2*(muAlpha-lambda1a)+kappa*(muBeta+lambda1b);
        numeratorBeta = phi1*(muBeta+lambda1b)+kappa*(muAlpha-lambda1a);
        if ((numeratorAlpha > 0) && (numeratorBeta < 0)) {
          alphaNew(0, j) = numeratorAlpha / denominator;
          betaNew(j, 0) = numeratorBeta / denominator;
        } else {
          numeratorAlpha = phi2*(muAlpha+lambda1a)+kappa*(muBeta-lambda1b);
          numeratorBeta = phi1*(muBeta-lambda1b)+kappa*(muAlpha+lambda1a);
          if ((numeratorAlpha < 0) && (numeratorBeta > 0)) {
            alphaNew(0, j) = numeratorAlpha / denominator;
            betaNew(j, 0) = numeratorBeta / denominator;
          } else {
            numeratorAlpha = phi2*(muAlpha+lambda1a)-kappa*(muBeta+lambda1b);
            numeratorBeta = phi1*(muBeta+lambda1b)-kappa*(muAlpha+lambda1a);
            if ((numeratorAlpha < 0) && (numeratorBeta < 0)) {
              alphaNew(0, j) = numeratorAlpha / denominator;
              betaNew(j, 0) = numeratorBeta / denominator;
            } else {
              numeratorAlpha = abs(muAlpha) - lambda1a;
              if ((numeratorAlpha > 0) && (phi1*abs(muBeta)-kappa*abs(muAlpha) <= phi1*lambda1b-kappa*lambda1a)) {
                alphaNew(0, j) = sgn(muAlpha) * numeratorAlpha / phi1;
                betaNew(j, 0) = 0.0;
              } else {
                numeratorBeta = abs(muBeta) - lambda1b;
                if ((numeratorBeta > 0) && (phi2*abs(muAlpha)-kappa*abs(muBeta) <= phi2*lambda1a-kappa*lambda1b)) {
                  alphaNew(0, j) = 0.0;
                  betaNew(j, 0) = sgn(muBeta) * numeratorBeta / phi2;
                } else {
                  alphaNew(0, j) = 0.0;
                  betaNew(j, 0) = 0.0;
                }
              }
            }
          }
        }
      }
    }
  }
  return std::make_tuple(alphaNew, betaNew);
};

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> upadteAlphaBetaPathwayNetwork(
    Eigen::MatrixXd laplacianMatrixA,
    Eigen::MatrixXd laplacianMatrixB,
    Eigen::MatrixXd alpha,
    Eigen::MatrixXd beta,
    Eigen::MatrixXd alphaStep1,
    Eigen::MatrixXd betaStep2,
    Eigen::MatrixXd tauAlpha,
    Eigen::MatrixXd tauBeta,
    double rho,
    double lambda1a,
    double lambda1b,
    double lambda2a,
    double lambda2b,
    double kappa,
    double lambda2aStar,
    double lambda2bStar
) {
  int p = alphaStep1.cols(), j;
  double Wa2 = lambda2a*lambda2aStar;
  double Wb2 = lambda2b*lambda2bStar;

  double phi1 = 2*kappa*Wa2+rho, phi2 = 2*kappa*Wb2+rho;

  double muAlpha, muBeta, denominator, numeratorAlpha, numeratorBeta;
  double Wa1, Wb1;

  Eigen::MatrixXd alphaNew = alpha, betaNew = beta;
  for (j = 0; j < p; ++j) {

    Wa1 = lambda2a*laplacianMatrixA.row(j).dot(alphaNew.row(0));
    Wb1 = lambda2b*laplacianMatrixB.row(j).dot(betaNew.col(0));

    muAlpha = -kappa*Wa1 + tauAlpha(0, j) + rho*alphaStep1(0, j);
    muBeta = -kappa*Wb1 + tauBeta(j, 0) + rho*betaStep2(j, 0);

    if (kappa == 0.0) {
      alphaNew(0, j) = softThreshold(muAlpha, lambda1a) / phi1;
      betaNew(j, 0) = softThreshold(muBeta, lambda1b) / phi2;
    } else {
      denominator = phi1*phi2-kappa*kappa;

      numeratorAlpha = phi2*(muAlpha-lambda1a)-kappa*(muBeta-lambda1b);
      numeratorBeta = phi1*(muBeta-lambda1b)-kappa*(muAlpha-lambda1a);
      if ((numeratorAlpha > 0) && (numeratorBeta > 0)) {
        alphaNew(0, j) = numeratorAlpha / denominator;
        betaNew(j, 0) = numeratorBeta / denominator;
      } else {
        numeratorAlpha = phi2*(muAlpha-lambda1a)+kappa*(muBeta+lambda1b);
        numeratorBeta = phi1*(muBeta+lambda1b)+kappa*(muAlpha-lambda1a);
        if ((numeratorAlpha > 0) && (numeratorBeta < 0)) {
          alphaNew(0, j) = numeratorAlpha / denominator;
          betaNew(j, 0) = numeratorBeta / denominator;
        } else {
          numeratorAlpha = phi2*(muAlpha+lambda1a)+kappa*(muBeta-lambda1b);
          numeratorBeta = phi1*(muBeta-lambda1b)+kappa*(muAlpha+lambda1a);
          if ((numeratorAlpha < 0) && (numeratorBeta > 0)) {
            alphaNew(0, j) = numeratorAlpha / denominator;
            betaNew(j, 0) = numeratorBeta / denominator;
          } else {
            numeratorAlpha = phi2*(muAlpha+lambda1a)-kappa*(muBeta+lambda1b);
            numeratorBeta = phi1*(muBeta+lambda1b)-kappa*(muAlpha+lambda1a);
            if ((numeratorAlpha < 0) && (numeratorBeta < 0)) {
              alphaNew(0, j) = numeratorAlpha / denominator;
              betaNew(j, 0) = numeratorBeta / denominator;
            } else {
              numeratorAlpha = abs(muAlpha) - lambda1a;
              if ((numeratorAlpha > 0) && (phi1*abs(muBeta)-kappa*abs(muAlpha) <= phi1*lambda1b-kappa*lambda1a)) {
                alphaNew(0, j) = sgn(muAlpha) * numeratorAlpha / phi1;
                betaNew(j, 0) = 0.0;
              } else {
                numeratorBeta = abs(muBeta) - lambda1b;
                if ((numeratorBeta > 0) && (phi2*abs(muAlpha)-kappa*abs(muBeta) <= phi2*lambda1a-kappa*lambda1b)) {
                  alphaNew(0, j) = 0.0;
                  betaNew(j, 0) = sgn(muBeta) * numeratorBeta / phi2;
                } else {
                  alphaNew(0, j) = 0.0;
                  betaNew(j, 0) = 0.0;
                }
              }
            }
          }
        }
      }
    }
  }
  return std::make_tuple(alphaNew, betaNew);
};

Eigen::MatrixXd updateGammaFunc(
    Eigen::MatrixXd betaStep2,
    Eigen::MatrixXd XtXInv,
    Eigen::MatrixXd XtM1,
    Eigen::MatrixXd XtY,
    double lambda1g
) {
  int p = XtXInv.cols(), j;
  Eigen::MatrixXd gammaTemp = XtY - XtM1 * betaStep2;
  for (j = 0; j < p; ++j) {
    gammaTemp(j, 0) = softThreshold(gammaTemp(j, 0), lambda1g);
  }
  return XtXInv * gammaTemp;
}

// [[Rcpp::export]]
Rcpp::List singleModalityAdmmFit(
    Eigen::Map<Eigen::MatrixXd> X,
    Eigen::Map<Eigen::MatrixXd> Y,
    Eigen::Map<Eigen::MatrixXd> M1,
    Eigen::Map<Eigen::MatrixXd> alphaInit,
    Eigen::Map<Eigen::MatrixXd> betaInit,
    Eigen::Map<Eigen::MatrixXd> gammaInit,
    double rho,
    double lambda1a,
    double lambda1b,
    double lambda1g,
    double lambda2a,
    double lambda2b,
    int penaltyType,
    Rcpp::List penaltyParameters,
    Eigen::Map<Eigen::MatrixXd> XtX,
    Eigen::Map<Eigen::MatrixXd> XtXInv,
    Eigen::Map<Eigen::MatrixXd> XtXPlusRhoInv,
    Eigen::Map<Eigen::MatrixXd> XtM1,
    Eigen::Map<Eigen::MatrixXd> M1tM1PlusRhoInv,
    Eigen::Map<Eigen::MatrixXd> M1tY,
    Eigen::Map<Eigen::MatrixXd> XtY,
    int maxIter,
    double tol,
    bool verbose,
    int verboseNumIter,
    int verboseNumAlpha,
    int verboseNumBeta,
    int verboseNumGamma
) {
  int p = M1.cols(), numColsX = X.cols();

  Eigen::MatrixXd laplacianMatrixA = Eigen::MatrixXd::Identity(1, 1);
  Eigen::MatrixXd laplacianMatrixB = Eigen::MatrixXd::Identity(1, 1);
  if (penaltyType == 2) {
    laplacianMatrixA = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixA"));
    laplacianMatrixB = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixB"));
  }

  double kappa = 0.0;
  if (penaltyType == 3) {
    kappa = Rcpp::as<double>(penaltyParameters("kappa"));
  }

  double lambda2aStar = 0.0, lambda2bStar = 0.0;
  if (penaltyType == 4) {
    kappa = Rcpp::as<double>(penaltyParameters("kappa"));
    lambda2aStar = Rcpp::as<double>(penaltyParameters("lambda2aStar"));
    lambda2bStar = Rcpp::as<double>(penaltyParameters("lambda2bStar"));
    laplacianMatrixA = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixA"));
    laplacianMatrixB = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixB"));
  }

  int iter = 0;
  bool converged = false;
  Eigen::MatrixXd alpha = alphaInit, beta = betaInit, gamma = gammaInit;
  Eigen::MatrixXd alphaStep1 = Eigen::MatrixXd::Zero(1, p), betaStep2 = Eigen::MatrixXd::Zero(p, 1);
  Eigen::MatrixXd tauAlpha = Eigen::MatrixXd::Zero(1, p), tauBeta = Eigen::MatrixXd::Zero(p, 1);
  Eigen::MatrixXd alphaStep1New, betaStep2New, alphaNew, betaNew, gammaNew, tauAlphaNew, tauBetaNew;
  while ((iter <= maxIter) && !converged) {
    ++iter;
    alphaStep1New = XtXPlusRhoInv * (XtM1 + rho*alpha - tauAlpha);
    betaStep2New = M1tM1PlusRhoInv * (M1tY - XtM1.transpose() * gamma + rho*beta - tauBeta);

    if (penaltyType == 1) {
      alphaNew = upadteAlphaElasticNet(alphaStep1New, tauAlpha, rho, lambda1a, lambda2a);
      betaNew = upadteBetaElasticNet(betaStep2New, tauBeta, rho, lambda1b, lambda2b);
    } else if (penaltyType == 2) {
      alphaNew = upadteAlphaNetwork(laplacianMatrixA, alpha, alphaStep1New, tauAlpha, rho, lambda1a, lambda2a);
      betaNew = upadteBetaNetwork(laplacianMatrixB, beta, betaStep2New, tauBeta, rho, lambda1b, lambda2b);
    } else if (penaltyType == 3) {
      std::tie(alphaNew, betaNew) = upadteAlphaBetaPathwayLasso(
        alphaStep1New, betaStep2New, tauAlpha, tauBeta,
        rho, lambda1a, lambda1b, kappa, lambda2a, lambda2b
      );
    }  else if (penaltyType == 4) {
      std::tie(alphaNew, betaNew) = upadteAlphaBetaPathwayNetwork(
        laplacianMatrixA, laplacianMatrixB,
        alpha, beta, alphaStep1New, betaStep2New, tauAlpha, tauBeta,
        rho, lambda1a, lambda1b, lambda2a, lambda2b, kappa, lambda2aStar, lambda2bStar
      );
    }

    gammaNew = updateGammaFunc(betaStep2New, XtXInv, XtM1, XtY, lambda1g);
    tauAlphaNew = tauAlpha + rho * (alphaStep1New - alphaNew);
    tauBetaNew = tauBeta + rho * (betaStep2New - betaNew);

    converged = ((gammaNew - gamma).array().pow(2).sum() < tol) &&
      ((alphaStep1New - alphaNew).array().pow(2).sum() < tol) &&
      ((alphaStep1New - alphaStep1).array().pow(2).sum() < tol) &&
      ((betaStep2New - betaNew).array().pow(2).sum() < tol) &&
      ((betaStep2New - betaStep2).array().pow(2).sum() < tol);

    if (verbose) {
      if (verbose && ((iter % verboseNumIter == 0) || converged)) {
        std::string isConvergedString = converged?"yes":"no";
        Rcpp::Rcout << std::fixed << std::setprecision(5);
        Rcpp::Rcout << "Iteration " << iter << ": is converged: " << isConvergedString;
        if (verboseNumGamma > 0) {
          printCoefficient(gammaNew.data(), "gamma", std::min(verboseNumGamma, numColsX));
        }
        if (verboseNumAlpha > 0) {
          printCoefficient(alphaNew.data(), "alpha", std::min(verboseNumAlpha, p));
        }
        if (verboseNumBeta > 0) {
          printCoefficient(betaNew.data(), "beta", std::min(verboseNumBeta, p));
        }
        Rcpp::Rcout << std::endl;
      }
    }

    alphaStep1 = alphaStep1New;
    betaStep2 = betaStep2New;
    alpha = alphaNew;
    beta = betaNew;
    gamma = gammaNew;
    tauAlpha = tauAlphaNew;
    tauBeta = tauBetaNew;
  }

  return Rcpp::List::create(
    Rcpp::Named("alpha") = alphaNew,
    Rcpp::Named("beta") = betaNew,
    Rcpp::Named("gamma") = gammaNew,
    Rcpp::Named("niter") = iter,
    Rcpp::Named("converged") = converged
  );
}
