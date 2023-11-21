// Copyright (C) 2022-2023             Ching-Chuan Chen
//
// This file is part of oneMKL.
//
// oneMKL.MatrixCal is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
//
// oneMKL.MatrixCal is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
#include <RcppEigen.h>

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

  if (Rf_ncols(X) != Rf_nrows(Y)) {
    Rcpp::stop("The number of rows of Y must be equal to the number of columns of X");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  if (Rf_ncols(Y) == 1) {
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

  if (Rf_nrows(X) != Rf_nrows(Y)) {
    Rcpp::stop("The number of rows of Y must be equal to the number of rows of X");
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
Eigen::MatrixXd fMatInv(SEXP X, bool is_sym_pd = false) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_ncols(X)) {
    Rcpp::stop("X must be a square matrix");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  if (is_sym_pd) {
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(XMtd.rows(), XMtd.cols());
    return XMtd.llt().solve(I);
  } else {
    return XMtd.inverse();
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd fMatChol(SEXP X){
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_ncols(X)) {
    Rcpp::stop("X must be a square matrix");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  return XMtd.llt().matrixU();
}
