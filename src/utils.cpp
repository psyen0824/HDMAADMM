#include <RcppEigen.h>
#include <Rcpp.h>

// [[Rcpp::export]]
SEXP cast_numeric(SEXP input) {
  if (!Rf_isReal(input)) {
    return Rf_coerceVector(input, REALSXP);
  } else {
    return input;
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd fMatProd(SEXP X, SEXP Y, bool is_X_symmetric = false) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  int *ydims;
  ydims = INTEGER(Rf_coerceVector(Rf_getAttrib(Y, R_DimSymbol), INTSXP));
  if (ydims[1] == 1) {
    Eigen::Map<Eigen::VectorXd> b = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * b;
    } else {
      return XMtd * b;
    }
  } else {
    Eigen::Map<Eigen::MatrixXd> Z = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * Z;
    } else {
      return XMtd * Z;
    }
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd fMatTransProd(SEXP X, SEXP Y, bool is_X_symmetric = false) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  int *ydims;
  ydims = INTEGER(Rf_coerceVector(Rf_getAttrib(Y, R_DimSymbol), INTSXP));
  if (ydims[1] == 1) {
    // to have the best performance
    Eigen::Map<Eigen::VectorXd> b = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * b;
    } else {
      return XMtd.transpose() * b;
    }
  } else {
    Eigen::Map<Eigen::MatrixXd> Z = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * Z;
    } else {
      return XMtd.transpose() * Z;
    }
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd fMatSolve(
    SEXP X, SEXP Y,
    bool is_sym_pd = false,
    bool is_invertible = false
) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
  if (is_sym_pd) {
    return XMtd.llt().solve(YMtd);
  } else if (is_invertible) {
    return XMtd.partialPivLu().solve(YMtd);
  } else {
    return XMtd.householderQr().solve(YMtd);
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd fMatInv(SEXP X, bool is_sym_pd = false) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  if (is_sym_pd) {
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(XMtd.rows(), XMtd.cols());
    return XMtd.llt().solve(I);
  } else {
    return XMtd.inverse();
  }
}
