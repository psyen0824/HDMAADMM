## Copyright (C) 2023        Ching-Chuan Chen, Pei-Shan Yen Chia-Wei Kuo
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

#' High-dimensional Single Modality Mediation Models
#'
#' @param X The matrix of independent variables (exposure/treatment/group).
#' @param Y The vector of dependent variable (outcome response).
#' @param M1 The single modality mediator.
#' @param rho The augmented Lagrangian parameter for ADMM.
#' @param lambda1g The L1-norm penalty for the direct effect. Default is \strong{10} to address overestimate issue.
#' @param lambda1a The L1-norm penalty for the effect between mediator and independent variables.
#' @param lambda1b The L1-norm penalty for the effect between mediator and dependent variable.
#' @param lambda2a The L2-norm penalty for the effect between mediator and independent variables.
#'   It's not used when Penalty is \code{PathwayLasso} which you can set it as 0.
#' @param lambda2b The L2-norm penalty for the effect between mediator and dependent variable.
#'   It's not used when Penalty is \code{PathwayLasso} which you can set it as 0.
#' @param penalty A string to specify the penalty. Default is \code{ElasticNet}. Possible methods are
#' Elastic Net (\code{ElasticNet}), Pathway Lasso (\code{PathywayLasso}), and  Network-constrained Penalty (\code{Network}).
#' @param penaltyParameterList
#' \itemize{
#'   \item Penalty=\code{ElasticNet} don't need other parameters.
#'   \item Penalty=\code{Network} needs two parameters.
#'   \itemize{
#'     \item laplacianMatrixA and laplacianMatrixB The L2-norm penalty for Network.
#'   }
#'   \item Penalty=\code{PathwayLasso} needs two parameters.
#'   \itemize{
#'     \item kappa The L1-norm penalty for Pathway Lasso.
#'     \item nu The L2-norm penalty for Pathway Lasso.
#'   }
#'   \item Penalty=\code{PathwayNetwork} needs five parameters.
#'   \itemize{
#'     \item kappa The L1-norm penalty for Pathway Network.
#'     \item lambda2aStar and lambda2bStar The L2-norm penalty for Pathway Network.
#'     \item laplacianMatrixA and laplacianMatrixB The L2-norm penalty for Pathway Network.
#'   }
#' }
#' @return A object, SingleModalityAdmm, with three elements.
#' \itemize{
#'   \item gamma: estimated direct effect.
#'   \item alpha: estimate effect between mediator and independent variables.
#'   \item beta : estimate effect between mediator and dependent variable.
#' }
#' @param SIS A logical value to specify whether to perform sure independence screening (SIS).
#' @param SISThreshold The threshold value for the target reduced dimension for mediators. The default is "2," which reduces the dimension to 2*n/log(n).
#' @param maxIter The maximum iterations. Default is \code{3000}.
#' @param tol The tolerance of convergence threshold. Default is \code{1e-3}.
#' @param verbose A logical value to specify whether to print the iteration process.
#' @param verboseOptions A list of values:
#' \itemize{
#'   \item \code{numIter}: The number of iterations to print.
#'   \item \code{numAlpha}: The number of \code{alpha} to print.
#'   \item \code{numBeta}: The number of \code{beta} to print.
#'   \item \code{numGamma}: The number of \code{gamma} to print.
#' }
#' @examples
#' ## Generate Empirical Data
#' simuData <- modalityMediationDataGen(seed = 20231201, generateLaplacianMatrix = TRUE)
#'
#' ## Parameter Estimation for ElasticNet penalty
#' modelElasticNet <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "ElasticNet"
#' )
#'
#' # fitted & predict
#' fitted(modelElasticNet)
#' predict(modelElasticNet, matrix(c(0, 1), ncol=1))
#'
#' ## Parameter Estimation for Pathway Lasso penalty
#' modelPathwayLasso <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "PathwayLasso", penaltyParameterList = list(kappa = 1, nu = 2)
#' )
#'
#' ## Parameter Estimation for Network penalty
#' modelNetwork <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "Network", penaltyParameterList = list(
#'     laplacianMatrixA = simuData$Info$laplacianMatrix, laplacianMatrixB = simuData$Info$laplacianMatrix
#'   )
#' )
#'
#' ## Parameter Estimation for Pathway Network penalty
#' modelPathwayNetwork <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "Network", penaltyParameterList = list(
#'     kappa = 1, lambda2aStar = 1, lambda2bStar = 1
#'     laplacianMatrixA = simuData$Info$laplacianMatrix, laplacianMatrixB = simuData$Info$laplacianMatrix
#'   )
#' )
#'
#' ## Parameter Estimation for Network penalty with a customized Laplacian matrix
#' set.seed(20231201)
#' p <- ncol(simuData$MediData$M1)
#' W <- matrix(0, nrow = p, ncol = p)
#' W[lower.tri(W)] <- runif(p*(p-1)/2, 0, 1)
#' W[upper.tri(W)] <- t(W)[upper.tri(W)]
#' diag(W) <- 1
#' L <- weightToLaplacian(W)
#' modelNetwork <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "Network", penaltyParameterList = list(laplacianMatrixA = L, laplacianMatrixB = L)
#' )
#'
#' ## With sure independence screening
#' ## Generate Empirical Data
#' simuData <- modalityMediationDataGen(n = 50, p = 1000, seed = 20231201)
#'
#' ## Parameter Estimation for ElasticNet penalty
#' modelElasticNetSIS <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "ElasticNet", SIS = TRUE
#' )
#' fitted(modelElasticNetSIS)
#' predict(modelElasticNetSIS, matrix(c(0, 1), ncol=1))
#' @importFrom stats lm coef cor
#' @export
singleModalityAdmm <- function(
    X, Y, M1,
    rho = 1, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b,
    penalty = "ElasticNet", penaltyParameterList = list(),
    SIS = FALSE, SISThreshold = 2,
    maxIter = 3000L, tol = 1e-3, verbose = FALSE,
    verboseOptions = list(numIter = 10L, numAlpha = 1L, numBeta = 1L, numGamma = 1L)
) {
  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X))
  }

  if (!is.matrix(Y)) {
    Y <- matrix(Y, nrow = length(Y))
  }

  if (ncol(X) != 1) {
    stop("The number of columns of X should be 1.")
  }

  if (ncol(Y) != 1) {
    stop("The number of columns of Y should be 1.")
  }

  if ((nrow(X) != length(Y)) || (nrow(X) != nrow(M1))) {
    stop("The length of Y should be equal to the rows of X and M1.")
  }

  if ((length(rho) > 1) || (length(lambda1a) > 1) || (length(lambda1b) > 1) || (length(lambda1g) > 1) || (length(lambda2a) > 1) || (length(lambda2b) > 1)) {
    stop("rho, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b should be scalar.")
  }

  if (any(is.na(X) | is.infinite(X))) {
    stop("X should be finite non-nan numeric matrix")
  }

  if (any(is.na(Y) | is.infinite(Y))) {
    stop("Y should be finite non-nan numeric matrix")
  }

  if (any(is.na(M1) | is.infinite(M1))) {
    stop("M1 should be finite non-nan numeric matrix")
  }

  if (is.na(rho) || is.infinite(rho)) {
    stop("rho should be finite non-nan numeric")
  }

  if (is.na(lambda1b) || is.infinite(lambda1b)) {
    stop("lambda1b should be finite non-nan numeric")
  }

  if (is.na(lambda1a) || is.infinite(lambda1a)) {
    stop("lambda1a should be finite non-nan numeric")
  }

  if (is.na(lambda2a) | is.infinite(lambda2a)) {
    stop("lambda2a should be finite non-nan numeric")
  }

  if (any(is.na(lambda2b) | is.infinite(lambda2b))) {
    stop("lambda2b should be finite non-nan numeric")
  }

  if (any(is.na(lambda1g) || is.infinite(lambda1g))) {
    stop("lambda1g should be finite non-nan numeric")
  }

  if (any(is.na(SISThreshold) || is.infinite(SISThreshold))) {
    stop("SISThreshold should be finite non-nan numeric")
  }

  if (any(is.na(maxIter) || is.infinite(maxIter))) {
    stop("maxIter should be finite non-nan numeric")
  }

  if (any(is.na(tol) || is.infinite(tol))) {
    stop("tol should be finite non-nan numeric")
  }

  defaultVerboseOptions <- list(numIter = 10L, numAlpha = 1L, numBeta = 1L, numGamma = 1L)
  for (nm in names(defaultVerboseOptions)) {
    if (nm %in% names(verboseOptions)) {
      verboseOptions[[nm]] <- as.integer(verboseOptions[[nm]])
    } else {
      verboseOptions[[nm]] <- defaultVerboseOptions[[nm]]
    }
  }

  for (nm in names(verboseOptions)) {
    if (any(is.na(verboseOptions[[nm]]) || is.infinite(verboseOptions[[nm]]))) {
      stop(sprintf("%s should be finite non-nan numeric", nm))
    }
  }

  if (penalty == "ElasticNet") {
    if (length(penaltyParameterList) >= 1) {
      warning("ElasticNet penalty does not take any parameters.")
    }
  } else if (penalty == "Network") {
    checkPenaltyParameterList(penaltyParameterList, c("laplacianMatrixA", "laplacianMatrixB"), penalty)
  } else if (penalty == "PathwayLasso") {
    checkPenaltyParameterList(penaltyParameterList, c("kappa", "nu"), penalty)
  } else if (penalty == "PathwayNetwork") {
    checkPenaltyParameterList(penaltyParameterList, c("kappa", "laplacianMatrixA", "laplacianMatrixB", "lambda2aStar", "lambda2bStar"), penalty)
  } else {
    stop("No such penalty.")
  }

  sisIndex <- 1:ncol(M1)
  if (SIS) {
    pearsonCors <- abs(cor(Y, M1))
    sisIndex <- which(rank(-pearsonCors) <= SISThreshold*nrow(M1)/log(nrow(M1)))
    if (penalty == "Network") {
      penaltyParameterList$laplacianMatrix <- penaltyParameterList$laplacianMatrix[sisIndex, sisIndex]
    }
  }

  YY <- scale(Y)
  Y.center <- attr(YY,"scaled:center")
  Y.scale <- attr(YY,"scaled:scale")
  XX <- scale(X)
  X.center <- attr(XX,"scaled:center")
  X.scale <- attr(XX,"scaled:scale")

  MM1 <- scale(M1[ , sisIndex])
  M1.center <- attr(MM1, "scaled:center")
  M1.scale <- attr(MM1, "scaled:scale")

  p <- ncol(MM1)

  gammaInit <- elasticNetFit(cbind(XX, MM1), YY, rep(0, p+1), lambda1g, 0., maxIter, tol)$coef[1]
  alphaInit <- matrix(sapply(1:p, function(i) {
    elasticNetFit(XX, MM1[ , i], rep(0, 1), lambda1a, lambda2a, maxIter, tol)$coef[1]
  }), nrow=1)
  betaInit <- matrix(sapply(1:p, function(i) {
    elasticNetFit(cbind(XX, MM1[ , i]), YY, rep(0, p+1), lambda1b, lambda2b, maxIter, tol)$coef[2]
  }), ncol=1)

  if (verbose) {
    if (verbose) {
      msg <- sprintf("Iteration 0: is converged: no, ")
      numGammaToPrint <- min(verboseOptions$numGamma, ncol(X))
      gammaMsg <- ""
      if (numGammaToPrint > 0) {
        gammaMsg <- paste0(sprintf("gamma[%i]: %.5f", 1:numGammaToPrint, gammaInit[1:numGammaToPrint]), collapse = ", ")
      }
      numAlphaToPrint <- min(verboseOptions$numAlpha, ncol(X))
      alphaMsg <- ""
      if (numAlphaToPrint > 0) {
        alphaMsg <- paste0(sprintf("alpha[%i]: %.5f", 1:numAlphaToPrint, alphaInit[1:numAlphaToPrint]), collapse = ", ")
      }
      numBetaToPrint <- min(verboseOptions$numBeta, ncol(X))
      betaMsg <- ""
      if (numBetaToPrint > 0) {
        betaMsg <- paste0(sprintf("beta[%i]: %.5f", 1:numBetaToPrint, betaInit[1:numBetaToPrint]), collapse = ", ")
      }
      cat(msg, gammaMsg, alphaMsg, betaMsg, "\n")
    }
  }

  XtX <- fMatTransProd(XX, XX)
  XtM1 <- fMatTransProd(XX, MM1)
  M1tM1PlusRhoInv <- fMatInv(fMatTransProd(MM1, MM1) + rho*diag(p))
  M1tY <- fMatTransProd(MM1, YY)
  XtY <- fMatTransProd(XX, YY)
  XtXInv <- fMatInv(XtX)
  XtXPlusRhoInv <- fMatInv(XtX + rho*diag(ncol(XX)), TRUE)

  penaltyType <- switch(penalty, ElasticNet = 1L, Network = 2L, PathwayLasso = 3L, PathwayNetwork = 4L)
  fitResult <- singleModalityAdmmFit(
    XX, YY, MM1, alphaInit, betaInit, gammaInit,
    rho, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b,
    penaltyType=penaltyType, penaltyParameters=penaltyParameterList,
    XtX=XtX, XtXInv=XtXInv, XtXPlusRhoInv=XtXPlusRhoInv, XtM1=XtM1,
    M1tM1PlusRhoInv=M1tM1PlusRhoInv, M1tY=M1tY,XtY=XtY,
    maxIter=as.integer(maxIter), tol=tol, verbose=verbose,
    verboseNumIter=verboseOptions$numIter, verboseNumAlpha=verboseOptions$numAlpha,
    verboseNumBeta=verboseOptions$numBeta, verboseNumGamma=verboseOptions$numGamma
  )

  if (fitResult$niter >= maxIter) {
    warning("Method does not converge!")
  }

  pTrue <- ncol(M1)
  alphaOut <- matrix(rep(0, pTrue), nrow=1)
  alphaOut[ , sisIndex] <- fitResult$alpha * M1.scale / X.scale
  betaOut <- matrix(rep(0, pTrue), ncol=1)
  betaOut[sisIndex, ] <- fitResult$beta * Y.scale / M1.scale
  gammaOut <- fitResult$gamma * Y.scale/X.scale

  out <- list(
    alpha = alphaOut,
    beta = betaOut,
    gamma = gammaOut,
    isConv = fitResult$converged,
    niter = fitResult$niter,
    interceptAlpha = colMeans(M1) - fMatProd(matrix(X.center, nrow=1), alphaOut),
    interceptBeta = as.numeric(
      Y.center - fMatProd(matrix(X.center, nrow=1), gammaOut) -
        fMatProd(fMatProd(matrix(X.center, nrow=1), alphaOut), betaOut)
    ),
    loglik = 0.0,
    fitted = matrix(rep(0, nrow(M1)), ncol=1)
  )

  out$loglik <- getLogLikelihood(
    X, Y, M1, out$alpha, out$beta, matrix(out$gamma, 1, 1))

  M1Hat <- sweep(fMatProd(X, out$alpha), 2, out$interceptAlpha, `+`)
  out$fitted <- out$interceptBeta + fMatProd(X, out$gamma) + fMatProd(M1Hat, out$beta)

  class(out) <- "SingleModalityAdmm"
  return(out)
}

#' Fitted Response of SingleModalityAdmm Fits
#'
#' @param object A fitted obejct of class inheriting from \code{SingleModalityAdmm}.
#' @param ... further arguments passed to or from other methods.
#' @method fitted SingleModalityAdmm
#' @return fitted.SingleModalityAdmm returns a vector which is fitted values.
#' @importFrom stats fitted
#' @export
fitted.SingleModalityAdmm <- function(object, ...) {
  return(object$fitted)
}

#' Predict Method for SingleModalityAdmm Fits
#'
#' @param object A fitted obejct of class inheriting from \code{SingleModalityAdmm}.
#' @param newdata Default is \code{NULL}. A matrix with variables to predict.
#' @param ... further arguments passed to or from other methods.
#' @method predict SingleModalityAdmm
#' @return predict.SingleModalityAdmm returns a vector which is the predicted values based on \code{newdata}.
#' @importFrom stats predict
#' @export
predict.SingleModalityAdmm <- function(object, newdata, ...) {
  if (is.null(newdata))
    return(fitted(object))
  M1Hat <- sweep(fMatProd(newdata, object$alpha), 2, object$interceptAlpha, `+`)
  return(object$interceptBeta + fMatProd(newdata, object$gamma) + fMatProd(M1Hat, object$beta))
}
