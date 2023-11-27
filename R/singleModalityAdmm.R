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

#' High-dimensional Single Modality Mediation Models
#'
#' @param X The matrix of independent variables (exposure/treatment/group).
#' @param Y The vector of dependent variable (outcome response).
#' @param M1 The single-modality mediator.
#' @param rho The augmented Lagrangian parameter for ADMM.
#' @param lambda1g The L1-norm penalty for the direct effect. Default is \strong{10} to adress overestimate issue.
#' @param lambda1a The L1-norm penalty for the effect between mediator and independent variables.
#' @param lambda1b The L1-norm penalty for the effect between mediator and dependent variable.
#' @param lambda2a The L2-norm penalty for the effect between mediator and independent variables.
#'   It's not used when Penalty is \code{PathwayLasso}.
#' @param lambda2b The L2-norm penalty for the effect between mediator and dependent variable.
#'   It's not used when Penalty is \code{PathwayLasso}.
#' @param penalty A string to specify the penalty. Default is \code{ElasticNet}. Possible methods are
#' Elastic Net (\code{ElasticNet}), Pathway Lasso (\code{PathywayLasso}), and  Network-constrained Penalty (\code{Network}).
#' @param penaltyParameterList
#' \itemize{
#'   \item Penalty=\code{PathwayLasso} needs two parameters.
#'   \itemize{
#'     \item kappa The L1-norm penalty for pathway Lasso.
#'     \item nu The L2-norm penalty for pathway Lasso.
#'   }
#'   \item Penalty=\code{Network} needs one parameter.
#'   \itemize{
#'     \item laplacianMatrix The Laplacian matrix applied on network penalty.
#'   }
#'   \item Penalty=\code{ElasticNet} don't need other parameters.
#' }
#' @return A object, SingleModalityAdmm, with three elements.
#' \itemize{
#'   \item gamma: estimated direct effect.
#'   \item alpha: estimate effect between mediator and independent variables.
#'   \item beta : estimate effect between mediator and dependent variables.
#' }
#' @param SIS A logical value to specify whether to perform sure independence screening (SIS).
#' @param SISThreshold The threshold value for the target reduced dimension for mediators. The default is "2," which reduces the dimension to 2*n/log(n).
#' @param maxIter The maximum iterations. Default is \code{3000}.
#' @param tol The tolerence of convergence threshold. Default is \code{1e-3}.
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
#'   penalty = "Network", penaltyParameterList = list(laplacianMatrix = simuData$Info$laplacianMatrix)
#' )
#'
#' ## Parameter Estimation for Network penalty with a customized Laplacian matrix
#' p <- 50
#' A <- matrix(rep(0, p*p), p, p)
#' A[1:10, 1:10] <- 1
#' A[1:10, 21:30] <- 1
#' A[21:30, 21:30] <- 1
#' A[21:30, 1:10] <- 1
#' A[11:20, 11:20] <- 1
#' A[11:20, 31:50] <- 1
#' A[31:50, 31:50] <- 1
#' A[31:50, 11:20] <- 1
#' diag(A) <- 0
#' L <- adjacencyToLaplacian(A)
#' modelNetwork <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
#'   penalty = "Network", penaltyParameterList = list(laplacianMatrix = L)
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
    X <- matrix(X, nrow = length(Y))
  }

  if (!is.matrix(Y)) {
    Y <- matrix(Y, nrow = length(Y))
  }

  defaultVerboseOptions <- list(numIter = 10L, numAlpha = 1L, numBeta = 1L, numGamma = 1L)
  for (nm in names(defaultVerboseOptions)) {
    if (nm %in% names(verboseOptions)) {
      verboseOptions[[nm]] <- as.integer(verboseOptions[[nm]])
    } else {
      verboseOptions[[nm]] <- defaultVerboseOptions[[nm]]
    }
  }

  sisIndex <- 1:ncol(M1)
  if (SIS) {
    pearsonCors <- abs(cor(Y, M1))
    sisIndex <- which(rank(-pearsonCors) <= SISThreshold*nrow(M1)/log(nrow(M1)))
  }

  YY <- scale(Y)
  Y.center <- attr(YY,"scaled:center")
  Y.scale <- attr(YY,"scaled:scale")
  XX <- scale(X)
  X.center <- attr(XX,"scaled:center")
  X.scale <- attr(XX,"scaled:scale")

  MM1 <- scale(M1[ , sisIndex])
  M1.center <- attr(MM1,"scaled:center")
  M1.scale <- attr(MM1,"scaled:scale")

  p <- ncol(MM1)

  gammaEst <- matrix(coef(lm(YY~XX+MM1))[2:(2+ncol(X)-1)], ncol=1)
  alphaEst <- coef(lm(MM1~XX))[2, , drop=FALSE]
  betaEst <- matrix(sapply(1:p, function(i) coef(lm(YY~XX+MM1[,i]))[3]), ncol=1)

  if (penalty == "Network") {
    if (!("laplacianMatrix" %in% names(penaltyParameterList))) {
      stop("penaltyParameterList should contains laplacianMatrix for Network penalty")
    }
  } else if (penalty == "PathwayLasso") {
    if (!("kappa" %in% names(penaltyParameterList)) || !("nu" %in% names(penaltyParameterList))) {
      stop("penaltyParameterList should contains kappa and nu for PathwayLasso penalty")
    }
  } else if (penalty == "ElasticNet") {
    # do nothing
  } else {
    stop("No such penalty.")
  }

  preCalcValues <- list(
    XtX = fMatTransProd(XX, XX),
    XtM1 = fMatTransProd(XX, MM1),
    M1tM1PlusRhoInv = fMatInv(fMatTransProd(MM1, MM1) + rho*diag(p)),
    M1tY = fMatTransProd(MM1, YY),
    XtY = fMatTransProd(XX, YY)
  )
  preCalcValues$XtXInv <- fMatInv(preCalcValues$XtX)
  preCalcValues$XtXPlusRhoInv <- fMatInv(preCalcValues$XtX + rho*diag(ncol(XX)), TRUE)

  penaltyType <- switch(penalty, ElasticNet = 1L, Network = 2L, PathwayLasso = 3L)
  fitResult <- singleModalityAdmmFit(
    XX, YY, MM1, alphaEst, betaEst, gammaEst,
    rho, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b,
    penaltyType, penaltyParameterList, as.integer(maxIter),
    XtX=preCalcValues$XtX, XtXInv = preCalcValues$XtXInv,
    XtXPlusRhoInv = preCalcValues$XtXPlusRhoInv, XtM1=preCalcValues$XtM1,
    M1tM1PlusRhoInv=preCalcValues$M1tM1PlusRhoInv, M1tY=preCalcValues$M1tY,
    XtY=preCalcValues$XtY, tol, verbose,
    verboseOptions$numIter, verboseOptions$numAlpha, verboseOptions$numBeta, verboseOptions$numGamma
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
    isConv = fitResult$niter < maxIter,
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
