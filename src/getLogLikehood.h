#include <RcppEigen.h>
#include "utility.h"

#ifndef HDMAADMM__LOGLIK__h
#define HDMAADMM__LOGLIK__h

Rcpp::List getLogLikelihood(
   Eigen::MatrixXd,
   Eigen::MatrixXd,
   Eigen::MatrixXd,
   Eigen::MatrixXd,
   Eigen::MatrixXd,
   Eigen::MatrixXd
);

double getObjective(
    Eigen::MatrixXd,
    Eigen::MatrixXd,
    Eigen::MatrixXd,
    Eigen::MatrixXd,
    Eigen::MatrixXd,
    Eigen::MatrixXd,
    int,
    double,
    double,
    double,
    double,
    double,
    double,
    Eigen::MatrixXd,
    Eigen::MatrixXd,
    double,
    double
);

#endif
