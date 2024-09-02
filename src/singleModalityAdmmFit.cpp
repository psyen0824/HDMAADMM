#include <optimization/LBFGSB.h>
#include <RcppEigen.h>
#include <tuple>
#include "utility.h"
#include "getLogLikehood.h"

// [[Rcpp::plugins(openmp)]]

class ElasticNetObjectiveGrad {
private:
  const Eigen::MatrixXd alphaStep1;
  const Eigen::MatrixXd betaStep2;
  const Eigen::MatrixXd tauAlpha;
  const Eigen::MatrixXd tauBeta;
  const int p;
  const double rho;
  const double lambda1a;
  const double lambda1b;
  const double lambda2a;
  const double lambda2b;
public:
  ElasticNetObjectiveGrad(
    const Eigen::MatrixXd alphaStep1_, const Eigen::MatrixXd betaStep2_,
    const Eigen::MatrixXd tauAlpha_, const Eigen::MatrixXd tauBeta_,
    const double rho_, const double lambda1a_, const double lambda1b_,
    const double lambda2a_, const double lambda2b_
  ): alphaStep1(alphaStep1_), betaStep2(betaStep2_), tauAlpha(tauAlpha_), tauBeta(tauBeta_), p(alphaStep1.cols()),
  rho(rho_), lambda1a(lambda1a_), lambda1b(lambda1b_), lambda2a(lambda2a_), lambda2b(lambda2b_)
  {}

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    Eigen::VectorXd alpha_pos(p), alpha_neg(p), beta_pos(p), beta_neg(p);
    for (int i = 0; i < p; ++i) {
      alpha_pos(i) = x(i);
      alpha_neg(i) = x(i + p);
      beta_pos(i) = x(i + 2*p);
      beta_neg(i) = x(i + 3*p);
    }

    Eigen::VectorXd alphaDiff = alphaStep1.row(0).transpose() - alpha_pos + alpha_neg;
    Eigen::VectorXd betaDiff = betaStep2.col(0) - beta_pos + beta_neg;

    double P2 = lambda1a * (alpha_pos.array().sum() + alpha_neg.array().sum()) + lambda1b * (beta_pos.array().sum() + beta_neg.array().sum());
    double P3 = lambda2a * (alpha_pos.dot(alpha_pos) + alpha_neg.dot(alpha_neg)) + lambda2b * (beta_pos.dot(beta_pos) + beta_neg.dot(beta_neg));
    double dual = tauAlpha.row(0).dot(alphaDiff) + tauBeta.col(0).dot(betaDiff);
    double f = P2 + P3 + dual + 0.5 * rho * (alphaDiff.dot(alphaDiff) + betaDiff.dot(betaDiff)); // rho/2 * ||x - z||^2

    Eigen::VectorXd gradP3_alpha(2 * p), gradP3_beta(2 * p);
    for (int i = 0; i < p; ++i) {
      gradP3_alpha(i) = 2.0 * lambda2a * alpha_pos(i);
      gradP3_alpha(i + p) = -2.0 * lambda2a * alpha_neg(i);
      gradP3_beta(i) = 2.0 * lambda2b * beta_pos(i);
      gradP3_beta(i + p) = -2.0 * lambda2b * beta_neg(i);
    }

    for (int i = 0; i < p; ++i) {
      grad(i) = lambda1a + gradP3_alpha(i) - tauAlpha(0, i) + rho * alphaDiff(i);
      grad(i + p) = lambda1a - gradP3_alpha(i+p) + tauAlpha(0, i) - rho * alphaDiff(i);
      grad(i + 2*p) = lambda1b + gradP3_alpha(i) - tauBeta(i, 0) + rho * betaDiff(i);
      grad(i + 3*p) = lambda1b - gradP3_beta(i+p) + tauBeta(i, 0) - rho * betaDiff(i);
    }
    return f;
  }
};


std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> updateAlphaBetaElasticNet(
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
    double lambda2b
) {
  int p = alpha.cols();

  Eigen::VectorXd x = Eigen::VectorXd::Zero(4*p);
  for (int i = 0; i < p; ++i) {
    if (alpha(0, i) >= 0) {
      x(i) = alpha(0, i);
    } else {
      x(i + p) = -alpha(0, i);
    }

    if (beta(i, 0) >= 0) {
      x(i + 2*p) = beta(i, 0);
    } else {
      x(i + 3*p) = -beta(i, 0);
    }
  }

  LBFGSpp::LBFGSBParam<double> param;
  param.epsilon        = 1e-5;
  param.epsilon_rel    = 1e-5;
  param.past           = 1;
  param.delta          = 1e-6;
  param.max_iterations = 10000;
  param.max_linesearch = 500;

  ElasticNetObjectiveGrad enog(alphaStep1, betaStep2, tauAlpha, tauBeta, rho, lambda1a, lambda1b, lambda2a, lambda2b);
  Eigen::VectorXd lb = Eigen::VectorXd::Zero(4*p);
  Eigen::VectorXd ub = Eigen::VectorXd::Zero(4*p);

  for (int i = 0; i < 4*p; ++i) {
    ub(i) = std::numeric_limits<double>::infinity();
  }


  LBFGSpp::LBFGSBSolver<double> solver(param);
  double fopt;
  int status = 0;
  try {
    solver.minimize(enog, x, fopt, lb, ub);
  } catch(const std::exception& e) {
    status = -1;
    Rcpp::warning(e.what());
  }

  Eigen::MatrixXd alphaNew(1, p), betaNew(p, 1);
  for (int i = 0; i < p; ++i) {
    alphaNew(0, i) = x(i) - x(i + p);
    betaNew(i, 0) = x(i + 2*p) - x(i + 3*p);
  }
  return std::make_tuple(alphaNew, betaNew);
}

Eigen::MatrixXd updateAlphaElasticNet(
    Eigen::MatrixXd alphaStep1,
    Eigen::MatrixXd tauAlpha,
    double rho,
    double lambda1a,
    double lambda2a
) {
  int p = alphaStep1.cols(), j;
  Eigen::MatrixXd alphaNew(1, p);

#if defined(_OPENMP)
#pragma omp for
#endif
  for (j = 0; j < p; ++j) {
    alphaNew(0, j) = softThreshold(tauAlpha(0, j) + rho*alphaStep1(0, j), lambda1a) / (lambda2a + rho);
  }
  return alphaNew;
};

Eigen::MatrixXd updateBetaElasticNet(
    Eigen::MatrixXd betaStep2,
    Eigen::MatrixXd tauBeta,
    double rho,
    double lambda1b,
    double lambda2b
) {
  int p = betaStep2.rows(), j;
  Eigen::MatrixXd betaNew(p, 1);

#if defined(_OPENMP)
#pragma omp for
#endif
  for (j = 0; j < p; ++j) {
    betaNew(j, 0) = softThreshold(tauBeta(j, 0) + rho*betaStep2(j, 0), lambda1b) / (lambda2b + rho);
  }
  return betaNew;
};

Eigen::MatrixXd updateAlphaNetwork(
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

Eigen::MatrixXd updateBetaNetwork(
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

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> updateAlphaBetaPathwayLasso(
    Eigen::MatrixXd alphaStep1,
    Eigen::MatrixXd betaStep2,
    Eigen::MatrixXd tauAlpha,
    Eigen::MatrixXd tauBeta,
    double rho,
    double lambda1a,
    double lambda1b,
    double lambda2a,
    double lambda2b,
    double kappa
) {
  int p = alphaStep1.cols(), j;
  // To do a fair comparison, we introduce lambda2a, lambda2b to represent L2 norm penalty
  double phi1 = 2*kappa*lambda2a+rho, phi2 = 2*kappa*lambda2b+rho;
  // Below is the old code: using nu
  // double phi1 = 2*kappa*nu+rho, phi2 = 2*kappa*nu+rho;
  double muAlpha, muBeta, denominator, numeratorAlpha, numeratorBeta;
  Eigen::MatrixXd alphaNew(1, p), betaNew(p, 1);

#if defined(_OPENMP)
#pragma omp for
#endif
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

/*
class PathwayNetworkObjectiveGrad: public Numer::MFuncGrad {
private:
  const Eigen::MatrixXd alphaStep1;
  const Eigen::MatrixXd betaStep2;
  const Eigen::MatrixXd tauAlpha;
  const Eigen::MatrixXd tauBeta;
  const Eigen::MatrixXd laplacianMatrixA;
  const Eigen::MatrixXd laplacianMatrixB;
  const int p;
  const double kappa;
  const double rho;
  const double lambda1a;
  const double lambda1b;
  const double lambda2a;
  const double lambda2b;
  const double lambda2aStar;
  const double lambda2bStar;
public:
  PathwayNetworkObjectiveGrad(
    const Eigen::MatrixXd alphaStep1_, const Eigen::MatrixXd betaStep2_,
    const Eigen::MatrixXd tauAlpha_, const Eigen::MatrixXd tauBeta_,
    const Eigen::MatrixXd laplacianMatrixA_, const Eigen::MatrixXd laplacianMatrixB_,
    const double kappa_, const double rho_, const double lambda1a_, const double lambda1b_,
    const double lambda2a_, const double lambda2b_, const double lambda2aStar_, const double lambda2bStar_
  ): alphaStep1(alphaStep1_), betaStep2(betaStep2_), tauAlpha(tauAlpha_), tauBeta(tauBeta_),
     laplacianMatrixA(laplacianMatrixA_), laplacianMatrixB(laplacianMatrixB_), p(alphaStep1.cols()),
     kappa(kappa_), rho(rho_), lambda1a(lambda1a_), lambda1b(lambda1b_), lambda2a(lambda2a_), lambda2b(lambda2b_),
     lambda2aStar(lambda2aStar_), lambda2bStar(lambda2bStar_)
  {}

  double f_grad(Numer::Constvec& x, Numer::Refvec grad) {
    Eigen::VectorXd alpha(p), beta(p);
    for (int i = 0; i < p; ++i) {
      alpha(i) = x(i);
      beta(i) = x(i + p);
    }

    Eigen::VectorXd alphaDiff = alphaStep1.row(0).transpose() - alpha;
    Eigen::VectorXd betaDiff = betaStep2.col(0) - beta;

    double P2 = lambda1a * alpha.array().abs().sum() + lambda1b * beta.array().abs().sum();
    double P3 = kappa * (
      (alpha.array() * beta.array()).abs().sum() +
        lambda2a * alpha.dot(laplacianMatrixA * alpha) +
        lambda2aStar * lambda2a * alpha.dot(alpha) +
        lambda2b * beta.dot(laplacianMatrixB * beta) +
        lambda2bStar * lambda2b * beta.dot(beta)
    );
    double dual = tauAlpha.row(0).dot(alphaDiff) + tauBeta.col(0).dot(betaDiff);
    const double f = P2 + P3 + dual + 0.5 * rho * (alphaDiff.dot(alphaDiff) + betaDiff.dot(betaDiff)); // rho/2 * ||x - z||^2

    Eigen::VectorXd gradP2_alpha = lambda1a * alpha.array().sign();
    Eigen::VectorXd gradP2_beta = lambda1b * beta.array().sign();
    Eigen::VectorXd gradP3_alpha = kappa * (beta.array().abs() * alpha.array().sign() +
      2.0 * lambda2a * (laplacianMatrixA * alpha).array() + 2.0 * lambda2a * lambda2aStar * alpha.array());
    Eigen::VectorXd gradP3_beta = kappa * (alpha.array().abs() * beta.array().sign() +
      2.0 * lambda2b * (laplacianMatrixB * beta).array() + 2.0 * lambda2b * lambda2bStar * beta.array());

    Eigen::VectorXd gradAlpha = gradP2_alpha + gradP3_alpha - tauAlpha.row(0).transpose() - rho * alphaDiff;
    Eigen::VectorXd gradBeta = gradP2_beta + gradP3_beta - tauBeta.col(0) - rho * betaDiff;
    for (int i = 0; i < p; ++i) {
      grad(i) = gradAlpha(i);
      grad(i + p) = gradBeta(i);
    }
    return f;
  }
};

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> updateAlphaBetaPathwayNetwork(
    Eigen::MatrixXd alpha,
    Eigen::MatrixXd beta,
    Eigen::MatrixXd laplacianMatrixA,
    Eigen::MatrixXd laplacianMatrixB,
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
  int p = alpha.cols();

  Eigen::VectorXd x(p + p);
  for (int i = 0; i < p; ++i) {
    x(i) = alpha(0, i);
    x(i + p) = beta(i, 0);
  }

  PathwayNetworkObjectiveGrad pnog(
    alphaStep1, betaStep2, tauAlpha, tauBeta, laplacianMatrixA, laplacianMatrixB,
    kappa, rho, lambda1a, lambda1b, lambda2a, lambda2b, lambda2aStar, lambda2bStar
  );

  double fopt;
  int status = Numer::optim_lbfgs(pnog, x, fopt, 5000, 1e-8, 1e-5);
  if (status < 0) {
    Rcpp::warning("algorithm did not converge");
  }

  Eigen::MatrixXd alphaNew(1, p), betaNew(p, 1);
  for (int i = 0; i < p; ++i) {
    alphaNew(0, i) = x(i);
    betaNew(i, 0) = x(i + p);
  }
  return std::make_tuple(alphaNew, betaNew);
};

*/

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

  // initialization for Network Penalty
  Eigen::MatrixXd laplacianMatrixA = Eigen::MatrixXd::Identity(1, 1);
  Eigen::MatrixXd laplacianMatrixB = Eigen::MatrixXd::Identity(1, 1);
  if (penaltyType == 2) {
    laplacianMatrixA = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixA"));
    laplacianMatrixB = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixB"));
  }

  // initialization for Pathway Lasso Penalty
  double kappa = 0.0;
  if (penaltyType == 3) {
    kappa = Rcpp::as<double>(penaltyParameters("kappa"));
  }

  // initialization for Pathway Network Penalty
  double lambda2aStar = 0.0, lambda2bStar = 0.0;
  if (penaltyType == 4) {
    kappa = Rcpp::as<double>(penaltyParameters("kappa"));
    lambda2aStar = Rcpp::as<double>(penaltyParameters("lambda2aStar"));
    lambda2bStar = Rcpp::as<double>(penaltyParameters("lambda2bStar"));
    laplacianMatrixA = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixA"));
    laplacianMatrixB = Rcpp::as<Eigen::MatrixXd>(penaltyParameters("laplacianMatrixB"));
  }

  // print objective for step 0 if verbose is true
  double obj = getObjective(
    X, Y, M1, alphaInit, betaInit, gammaInit,
    penaltyType, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b,
    kappa, laplacianMatrixA, laplacianMatrixB, lambda2aStar, lambda2bStar
  );
  if (verbose) {
    Rcpp::Rcout << std::fixed << std::setprecision(5);
    Rcpp::Rcout << "Iteration 0: is converged: no; Objective: " << obj;
    Rcpp::Rcout << "; gammaConv: 0, alphaConv: 0, betaConv: 0" << std::endl;
  }

  int iter = 0;
  bool gammaConv, alphaConv, betaConv, converged = false;
  Eigen::MatrixXd alpha = alphaInit, beta = betaInit, gamma = gammaInit;
  Eigen::MatrixXd alphaStep1 = Eigen::MatrixXd::Zero(1, p), betaStep2 = Eigen::MatrixXd::Zero(p, 1);
  Eigen::MatrixXd tauAlpha = Eigen::MatrixXd::Zero(1, p), tauBeta = Eigen::MatrixXd::Zero(p, 1);
  Eigen::MatrixXd alphaStep1New, betaStep2New, alphaNew, betaNew, gammaNew, tauAlphaNew, tauBetaNew;
  while ((iter <= maxIter) && !converged) {
    ++iter;
    // ADMM step 1: likelihood part
    alphaStep1New = XtXPlusRhoInv * (XtM1 + rho*alpha - tauAlpha);
    betaStep2New = M1tM1PlusRhoInv * (M1tY - XtM1.transpose() * gamma + rho*beta - tauBeta);
    gammaNew = updateGammaFunc(betaStep2New, XtXInv, XtM1, XtY, lambda1g);

    // ADMM step 2: penalty term part
    if (penaltyType == 1) {
      // alphaNew = updateAlphaElasticNet(alphaStep1New, tauAlpha, rho, lambda1a, lambda2a);
      // betaNew = updateBetaElasticNet(betaStep2New, tauBeta, rho, lambda1b, lambda2b);
      std::tie(alphaNew, betaNew) = updateAlphaBetaElasticNet(
        alpha, beta, alphaStep1New, betaStep2New, tauAlpha, tauBeta, rho, lambda1a, lambda1b, lambda2a, lambda2b
      );
    } else if (penaltyType == 2) {
      alphaNew = updateAlphaNetwork(laplacianMatrixA, alpha, alphaStep1New, tauAlpha, rho, lambda1a, lambda2a);
      betaNew = updateBetaNetwork(laplacianMatrixB, beta, betaStep2New, tauBeta, rho, lambda1b, lambda2b);
    } else if (penaltyType == 3) {
      std::tie(alphaNew, betaNew) = updateAlphaBetaPathwayLasso(
        alphaStep1New, betaStep2New, tauAlpha, tauBeta,
        rho, lambda1a, lambda1b, lambda2a, lambda2b, kappa
      );
    }  else if (penaltyType == 4) {
      /*
      std::tie(alphaNew, betaNew) = updateAlphaBetaPathwayNetwork(
        alpha, beta, laplacianMatrixA, laplacianMatrixB, alphaStep1New, betaStep2New, tauAlpha, tauBeta,
        rho, lambda1a, lambda1b, lambda2a, lambda2b, kappa, lambda2aStar, lambda2bStar
      );
      */
    }

    // ADMM step 3: dual update
    tauAlphaNew = tauAlpha + rho * (alphaStep1New - alphaNew);
    tauBetaNew = tauBeta + rho * (betaStep2New - betaNew);

    // check convergence
    gammaConv = ((gammaNew - gamma).array().pow(2).sum() < tol);
    alphaConv = ((alphaStep1New - alphaNew).array().pow(2).sum() < tol) &&
      ((alphaStep1New - alphaStep1).array().pow(2).sum() < tol);
    betaConv = ((betaStep2New - betaNew).array().pow(2).sum() < tol) &&
      ((betaStep2New - betaStep2).array().pow(2).sum() < tol);
    converged = (gammaConv && alphaConv && betaConv);

    // print coefs and objective if verbose is true
    if (verbose) {
      if (verbose && ((iter % verboseNumIter == 0) || converged)) {
        // calculate objective value
        obj = getObjective(
          X, Y, M1, alphaNew, betaNew, gammaNew,
          penaltyType, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b,
          kappa, laplacianMatrixA, laplacianMatrixB, lambda2aStar, lambda2bStar
        );
        // convert flag to string
        std::string isConvergedString = converged?"yes":"no";

        // print debugging information
        Rcpp::Rcout << "Iteration " << iter << ": is converged: " << isConvergedString << "; Objective: " << obj;
        Rcpp::Rcout << "; gammaConv: " << gammaConv << ", alphaConv: " << alphaConv << ", betaConv: " << betaConv;
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

    // update variables
    alphaStep1 = alphaStep1New;
    betaStep2 = betaStep2New;
    alpha = alphaNew;
    beta = betaNew;
    gamma = gammaNew;
    tauAlpha = tauAlphaNew;
    tauBeta = tauBetaNew;
  }

  // calculate log likelihood
  Rcpp::List logLikelihoodList = getLogLikelihood(X, Y, M1, alphaNew, betaNew, gammaNew);

  return Rcpp::List::create(
    Rcpp::Named("alpha") = alphaNew,
    Rcpp::Named("beta") = betaNew,
    Rcpp::Named("gamma") = gammaNew,
    Rcpp::Named("logLik") = logLikelihoodList,
    Rcpp::Named("niter") = iter,
    Rcpp::Named("converged") = converged
  );
}
