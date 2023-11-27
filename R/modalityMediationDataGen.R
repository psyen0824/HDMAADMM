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

#' @importFrom dqrng dqrnorm
rmvnorm <- function (n, mean, sigma = diag(length(mean))) {
  return(sweep(
    fMatProd(matrix(dqrnorm(n * ncol(sigma)), nrow = n), fMatChol(sigma)), 2, mean, "+"
  ))
}

#' Data Generation for High-Dimensional Mediation Model
#'
#' @param n The number of subjects for the high-dimensional mediation model)
#' @param p The number of high-dimensional mediators.
#' @param sigmaY The argument "sigmaY" represents the standard deviation (SD) of the error distribution for the dependent variable.
#' @param sizeNonZero The number of nonzero mediators. Here, we provide simulated scenarios that could produce large, medium,
#'  and small mediated effects, generating from a normal distribution.
#' @param alphaMean,alphaSd The mean and SD vector of the effect between the mediator and independent variable.
#' @param betaMean,betaSd The mean and SD vector of the effect between the mediator and dependent variable.
#' @param sigmaM1 The covariance matrix of the error distribution among mediators. Default is \code{diag(p)}.
#' @param gamma The true value of direct effect.
#' @param generateLaplacianMatrix A logical value to specify whether to generate Laplacian matrix for network penalty.
#' @param seed The random seed. Default is NULL to use the current seed.
#' @return A object with three elements.
#' \itemize{
#'   \item MediData: The simulated data for high-dimensional mediation model.
#'   \item MediPara: The true value for mediated effect and direct effect.
#'   \item Info : The output includes random seed, parameter setting, and Laplacian matrix for generating mediation model.
#' }
#' @examples
#' simuData <- modalityMediationDataGen(seed = 20231201)
#' @importFrom dqrng dqrunif dqrnorm dqset.seed
#' @export
modalityMediationDataGen <- function(
    n = 100,
    p = 50,
    sigmaY = 1,
    sizeNonZero = c(3, 3, 4),
    alphaMean = c(6, 4, 2),
    alphaSd = 0.1,
    betaMean = c(6, 4, 2),
    betaSd = 0.1,
    sigmaM1 = NULL,
    gamma = 3,
    generateLaplacianMatrix = FALSE,
    seed = 20231201
) {
  if (is.null(sigmaM1)) {
    sigmaM1 <- diag(p)
  }

  if ((ncol(sigmaM1) != nrow(sigmaM1)) || (ncol(sigmaM1) != p) || !isSymmetric.matrix(sigmaM1)) {
    stop("sigmaM1 must be a symmetric square matrix which is also positive definitie")
  }

  if (length(alphaSd) == 1) {
    alphaSd <- rep(alphaSd, length(alphaMean))
  }

  if (length(betaSd) == 1) {
    betaSd <- rep(betaSd, length(betaMean))
  }

  sizeVec <- c(length(sizeNonZero), length(alphaMean), length(alphaSd), length(betaMean), length(betaSd))
  if (length(unique(sizeVec)) > 1) {
    stop("sizeNonZero, alphaMean, alphaSd, betaMean, and betaSd must have the same length")
  }

  # random seed to generate multiple datasets
  dqset.seed(as.integer(seed))

  # parameters for causal effect
  # direct effect
  gamma <- matrix(gamma, ncol = 1)

  # mediated effect
  alpha <- matrix(
    c(
      unlist(mapply(
        function(n, m, s) dqrnorm(n, m, s),
        sizeNonZero, alphaMean, alphaSd,
        SIMPLIFY = FALSE
      )),
      unlist(mapply(
        function(n, m, s) dqrnorm(n, m, s),
        sizeNonZero, alphaMean, alphaSd,
        SIMPLIFY = FALSE
      )),
      rep(0, sum(sizeNonZero)),
      rep(0, p - 3*sum(sizeNonZero))
    ),
    nrow = 1
  )
  beta <- matrix(
    c(
      unlist(mapply(
        function(n, m, s) dqrnorm(n, m, s),
        sizeNonZero, betaMean, betaSd,
        SIMPLIFY = FALSE
      )),
      rep(0, sum(sizeNonZero)),
      unlist(mapply(
        function(n, m, s) dqrnorm(n, m, s),
        sizeNonZero, betaMean, betaSd,
        SIMPLIFY = FALSE
      )),
      rep(0, p - 3*sum(sizeNonZero))
    ),
    ncol = 1
  )

  # X = Binary Exposure/Treatment/group: generate X from a Bernoulli distribution
  X <- matrix(as.numeric(dqrunif(n) > 0.5), nrow = n, byrow = TRUE)

  # M1 = Mediator ; number of mediators = p1
  M1 <- fMatProd(X, alpha) + rmvnorm(n = n, mean = rep(0, p), sigma = sigmaM1)

  # Y = Continuous Outcome Response
  Y <- fMatProd(X, gamma) + fMatProd(M1, beta) + dqrnorm(n = n, mean = 0, sd = sigmaY)

  out <- list(
      MediData = list(X = X, M1= M1,  Y=Y),
      MediPara = list(alpha = alpha, beta = beta, gamma = gamma),
      Info = list(
        parameters = list(
          sigmaY = sigmaY,
          sizeNonZero = sizeNonZero,
          alphaMean = alphaMean,
          alphaSd = alphaSd,
          betaMean = betaMean,
          betaSd = betaSd,
          sigmaM1 = sigmaM1
        ),
        trueValue = list(gamma = gamma),
        laplacianMatrix = NULL,
        seed = seed
      )
    )

  if (generateLaplacianMatrix) {
    A <- matrix(0, p, p)
    for (i in 1:p) {
      A[i, i] <- summary(lm(Y ~ M1[ ,i]))$r.squared
    }
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        A[i, j] <- summary(lm(Y ~ M1[ ,i] + M1[ ,j]))$r.squared
      }
    }
    A[lower.tri(A)] <- t(A)[lower.tri(A)]
    out$Info$laplacianMatrix <- adjacencyToLaplacian(A)
  }
  return(out)
}
