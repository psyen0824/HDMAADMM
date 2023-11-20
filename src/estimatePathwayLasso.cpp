// Copyright (C) 2022-2023             Ching-Chuan Chen, Pei-Shan Yen
//
// This file is part of HDMAAMDD.
//
// HDMAAMDD is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
//
// HDMAAMDD is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
#include <RcppEigen.h>
#include "utility.h"

// [[Rcpp::export]]
Rcpp::List estimatePathwayLasso(
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
    double kappa,
    double nu
) {
  int p = M1.cols(), p2 = X.cols(), j;
  Eigen::MatrixXd alphaStep1 = (XtX + Eigen::MatrixXd::Identity(p2, p2)).llt().solve(XtM1 + rho*alpha - tauAlpha);
  Eigen::MatrixXd betaStep2 = M1tM1PlusRhoInv * (M1tY - XtM1.transpose() * gamma + rho*beta - tauBeta);

  double phi1 = 2*kappa*nu+rho, phi2 = 2*kappa*nu+rho;
  double muAlpha, muBeta, denominator, numeratorAlpha, numeratorBeta;
  Eigen::MatrixXd alphaNew = alpha, betaNew = beta;
  for (j = 0; j < p; j++) {
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
