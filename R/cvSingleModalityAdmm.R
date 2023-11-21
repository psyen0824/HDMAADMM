## Copyright (C) 2023        Ching-Chuan Chen, Pei-Shan Yen
##
## This file is part of HDMAADMM.
##
## HDMAADMM is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## HDMAADMM is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

#' Cross Validation for High-dimensional Single Mediation Models
#'
#' @param X The matrix of independent variables (exposure/treatment/group).
#' @param Y The vector of dependent variable (outcome response).
#' @param M1 The single-modality mediator.
#' @param nfolds The number of folds. The default is 10. Although nfolds can be as large as the sample size (leave-one-out CV),
#'  it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' @param typeMeasure Default is "rmse".
#' @param rho,lambda1g,lambda1a,lambda1b,lambda2a,lambda2b,penaltyParameterList Allow to put sequences for each parameter. Please refer to the function, \code{\link{singleModalityAdmm}} for the details.
#' @param penalty,SIS,SISThreshold,maxIter,tol,verbose,debug Please refer to the function, \code{\link{singleModalityAdmm}}.
#' @examples
#' ## Generate Empirical Data
#' simuData <- modalityMediationDataGen(seed = 20231201)
#'
#' ## Parameter Estimation for ElasticNet penalty
#' modelElasticNet <- cvSingleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   numFolds = 5, typeMeasure = "rmse",
#'   rho = 1, lambda1a = c(0.1, 0.5, 1), lambda1b = c(0.1, 0.3),
#'   lambda1g = 2, lambda2a = c(0.5, 1), lambda2b = c(0.5, 1),
#'   penalty = "ElasticNet"
#' )
#'
#' ## Parameter Estimation for Pathway Lasso penalty
#' modelPathwayLasso <- cvSingleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   numFolds = 5, typeMeasure = "rmse",
#'   rho = 1, lambda1a = c(0.1, 0.5, 1), lambda1b = c(0.1, 0.3),
#'   lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "PathwayLasso", penaltyParameterList = list(kappa = 1, nu = 2)
#' )
#' @export
cvSingleModalityAdmm <- function(
    X, Y, M1, numFolds = 10, typeMeasure = "rmse",
    lambda1a, lambda1b, lambda1g, lambda2a, lambda2b, rho = 1,
    penalty = "ElasticNet", penaltyParameterList = list(),
    SIS = FALSE, SISThreshold = 2,
    maxIter=3000, tol=1e-4, verbose = FALSE, debug = FALSE
) {
  if (numFolds < 3) {
    stop("numFolds must be greater than 3; suggests numFolds = 10.")
  }

  n <- nrow(M1)
  cvVec <- rep(1:numFolds, ceiling(n / numFolds), length.out = n)
  measureFunc <- function(Y, yHat, typeMeasure) {
    switch(
      typeMeasure,
      rmse = sqrt(mean((Y - yHat)^2))
    )
  }
  if (penalty == "PathwayLasso") {
    measureResultMatrix <- as.matrix(
      expand.grid(
        rho=rho, lambda1a=lambda1a, lambda1b=lambda1b, lambda1g=lambda1g,
        kappa=penaltyParameterList$kappa, nu=penaltyParameterList$nu,
        measure = 0.0
      )
    )

    for (j in 1:nrow(measureResultMatrix)) {
      foldResults <- vector("numeric", numFolds)
      for (f in 1:numFolds) {
        cvModel <- singleModalityAdmm(
          X = X[cvVec != f, , drop = FALSE], Y = Y[cvVec != f, , drop = FALSE], M1 = M1[cvVec != f, , drop = FALSE],
          rho = measureResultMatrix[j, 1], lambda1a = measureResultMatrix[j, 2],
          lambda1b = measureResultMatrix[j, 3], lambda1g = measureResultMatrix[j, 4],
          lambda2a = 1, lambda2b = 1, penalty = "PathwayLasso",
          penaltyParameterList = list(kappa = measureResultMatrix[j, 5], nu = measureResultMatrix[j, 6])
        )
        foldResults[f] = measureFunc(Y[cvVec == f, , drop = FALSE], predict(cvModel, X[cvVec == f, , drop = FALSE]), typeMeasure)
      }
      measureResultMatrix[j, 7] = mean(foldResults)
    }
  } else {
    measureResultMatrix <- as.matrix(
      expand.grid(
        rho=rho, lambda1a=lambda1a, lambda1b=lambda1b,
        lambda1g=lambda1g, lambda2a=lambda2a, lambda2b=lambda2b,
        measure = 0.0
      )
    )

    for (j in 1:nrow(measureResultMatrix)) {
      foldResults <- vector("numeric", numFolds)
      for (f in 1:numFolds) {
        cvModel <- singleModalityAdmm(
          X = X[cvVec != f, , drop = FALSE], Y = Y[cvVec != f, , drop = FALSE], M1 = M1[cvVec != f, , drop = FALSE],
          rho = measureResultMatrix[j, 1], lambda1a = measureResultMatrix[j, 2],
          lambda1b = measureResultMatrix[j, 3], lambda1g = measureResultMatrix[j, 4],
          lambda2a = measureResultMatrix[j, 5], lambda2b = measureResultMatrix[j, 6],
          penalty = penalty, penaltyParameterList = penaltyParameterList
        )
        foldResults[f] = measureFunc(Y[cvVec == f, , drop = FALSE], predict(cvModel, X[cvVec == f, , drop = FALSE]), typeMeasure)
      }
      measureResultMatrix[j, 7] = mean(foldResults)
    }
  }

  class(measureResultMatrix) <- "cvSingleModalityAdmm"
  return(measureResultMatrix)
}
