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

rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean))) {
  R <- chol(sigma, pivot = TRUE)
  R[, order(attr(R, "pivot"))]

  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% R
  retval <- sweep(retval, 2, mean, "+")
  retval
}

#' Data Generation for High-Dimensional Mediation Model
#'
#' @param n The number of subjects for the high-dimensional mediation model)
#' @param p The number of high-dimensional mediators.
#' @param parameters The true value of high-dimensional mediator model:
#'  The argument "SigmaY" represents the standard deviation (SD) of the error distribution for the dependent variable.
#'  The argument "SigmaM1" represents the covariance matrix of the error distribution among mediators.
#'  The argument "sizeNonZero" represents the number of nonzero mediators.
#'  Here, we provide simulated scenarios that could produce large, medium,
#'  and small mediated effects, generating from a normal distribution.
#'  The argument "alphaMean" and "alphaSD" represent the mean and SD vector of
#'  the effect between the mediator and independent variable.
#'  The argument "betaMean" and "betaSD" represent the mean and SD vector of
#'  the effect between the mediator and dependent variable.
#'  The argument "gamma" represents the true value of direct effect.
#' @param seed The random seed. Default is NULL to use the current seed.
#' @return A gen_data object with three elements.
#' \itemize{
#'   \item Medi_Data: The simulated data for high-dimensional mediation model.
#'   \item Medi_Para: The true value for mediated effect and direct effect.
#'   \item Info : The output includes random seed and parameter setting for generating mediation model.
#' }
#' @examples
#' simuData <- modalityMediationDataGen(seed = 20231201)
#' @export
modalityMediationDataGen <- function(
  n = 100, p = 50,
  parameters = list(
    sigmaY = 1,
    sigmaM1 = diag(p),
    sizeNonZero = c(3, 3, 4),
    alphaMean = c(6, 4, 2),
    alphaSd = c(0.1, 0.1, 0.1),
    betaMean = c(6,4,2),
    betaSd = c(0.1, 0.1, 0.1)
  ),
  trueValue = list(gamma = 3),
  seed = 20231117
  ) {
  # random seed to generate multiple datasets
  set.seed(seed)

  # parameters for causal effect
  # direct effect
  gamma <- matrix(trueValue$gamma, ncol = 1)

  # mediated effect
  alpha <- matrix(
    c(
      unlist(mapply(
        function(n, m, s) rnorm(n, m, s),
        parameters$sizeNonZero, parameters$alphaMean, parameters$alphaSd,
        SIMPLIFY = FALSE
      )),
      unlist(mapply(
        function(n, m, s) rnorm(n, m, s),
        parameters$sizeNonZero, parameters$alphaMean, parameters$alphaSd,
        SIMPLIFY = FALSE
      )),
      rep(0, sum(parameters$sizeNonZero)),
      rep(0, p - 3*sum(parameters$sizeNonZero))
    ),
    nrow = 1
  )
  beta <- matrix(
    c(
      unlist(mapply(
        function(n, m, s) rnorm(n, m, s),
        parameters$sizeNonZero, parameters$betaMean, parameters$betaSd,
        SIMPLIFY = FALSE
      )),
      rep(0, sum(parameters$sizeNonZero)),
      unlist(mapply(
        function(n, m, s) rnorm(n, m, s),
        parameters$sizeNonZero, parameters$betaMean, parameters$betaSd,
        SIMPLIFY = FALSE
      )),
      rep(0, p - 3*sum(parameters$sizeNonZero))
    ),
    ncol = 1
  )

  # X = Binary Exposure/Treatment/group: generate X from a Bernoulli distribution
  X <- as.matrix(rbinom(n = n, size = 1, prob = 0.5), nrow = n, byrow = TRUE)

  # M1 = Structural Mediator ; number of mediators = p1
  M1 <- X %*% alpha + rmvnorm(n = n, mean = rep(0, p), sigma = parameters$sigmaM1)

  # Y = Continuous Outcome Response
  Y <- X %*% gamma + M1 %*% tau + rnorm(n = n, mean = 0, sd = parameters$sigmaY)

  return(
    list(
      MediData = list(X = X, M1= M1,  Y=Y),
      MediPara = list(alpha = alpha, beta = beta, gamma = gamma),
      Info = list(parameters = parameters, trueValue = trueValue, seed = seed)
    )
  )
}
