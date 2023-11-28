#include <RcppEigen.h>
#include "utility.h"

// [[Rcpp::export]]
Rcpp::List elasticNetFit(
    Eigen::Map<Eigen::MatrixXd> X,
    Eigen::Map<Eigen::VectorXd> y,
    Eigen::Map<Eigen::VectorXd> coefInit,
    double lambda1,
    double lambda2,
    int maxIter = 3000,
    double tol = 1e-3,
    bool verbose = false,
    int verboseNumIter = 10,
    int verboseNumCoef = 1
) {
  int n = X.rows(), p = X.cols(), iter = 0, j;
  bool converged = false;
  double obj = 0.0, objNew = 9999.0, l2p1 = lambda2 + 1;
  Eigen::VectorXd rr = y - X * coefInit;
  Eigen::VectorXd coef = coefInit, coefNew = coefInit;
  while ((iter <= maxIter) && !converged) {
    iter += 1;
    obj = objNew;
    coef = coefNew;
    for (j = 0; j < p; j++) {
      coefNew(j) = softThreshold(X.col(j).dot(rr) / n + coef(j) * l2p1, lambda1) / l2p1;
      rr -= (coefNew(j) - coef(j)) * X.col(j);
    }

    objNew = rr.array().pow(2).sum()/n + lambda1 * coefNew.array().abs().sum() + lambda2 * coefNew.dot(coefNew) / 2.0;
    converged = (std::abs(objNew - obj) / objNew < tol);

    if (verbose) {
      if (verbose && ((iter % verboseNumIter == 0) || converged)) {
        std::string isConvergedString = converged?"yes":"no";
        Rcpp::Rcout << std::fixed << std::setprecision(5);
        Rcpp::Rcout << "Iteration " << iter << ": is converged: " << isConvergedString;
        if (verboseNumCoef > 0) {
          printCoefficient(coefNew.data(), "coef", std::min(verboseNumCoef, p));
        }
        Rcpp::Rcout << std::endl;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("coef") = coefNew,
    Rcpp::Named("niter") = iter,
    Rcpp::Named("converged") = converged
  );
}
