#' High-dimensional Mediation Models with Pathway Lasso Penalty.
#'
#' @param X The matrix of independent variables (exposure/treatment/group).
#' @param Y The vector of dependent variable (outcome response).
#' @param M1 The single-modality mediator.
#' @param alpha The initial value for the effect between mediator and independent variables.
#' @param beta The initial value for the effect between mediator and dependent variable.
#' @param gamma The initial value for the direct effect.
#' @param tauAlpha The iteration parameter for ADMM - for the effect between mediator and independent variable.
#' @param tauBeta The iteration parameter for ADMM - for the effect between mediator and depenent variable.
#' @param rho The augmented Lagrangian parameter for ADMM.
#' @param lambda1g The L1-norm penalty for the direct effect. Default is "10" to adress overestimate issue.
#' @param lambda1a The L1-norm penalty for the effect between mediator and independent variables.
#' @param lambda1b The L1-norm penalty for the effect between mediator and dependent variable.
#' @param lambda2a The L2-norm penalty for the effect between mediator and independent variables.
#' @param lambda2b The L2-norm penalty for the effect between mediator and dependent variable.
#' @param penalty A string to specify the penalty. Default is "Lasso". Possible methods are
#' Lasso ("Lasso"), Elastic Net ("ElasticNet"), Pathway Lasso ("PathwayLasso"), and  Network-constrained Penalty ("Network").
#' @param SIS A logical value to specify whether to perform sure independence screening (SIS).
#' @param SISThreshold The threshold value for the target reduced dimension for mediators. The default is "2," which reduces the dimension to 2*n/log(n).
#' @param XtX Input the multiplication of matrices `X^(t)` and `X`, i.e., `X^(t)X`.
#' @param XtM1 Input the multiplication of matrices `X^(t)` and `M1`, i.e., `X^(t)M1`.
#' @param M1tM1PlusRhoInv Input the matrix to estimate the effect between mediator and independned variables.
#' @param M1tY Input the multiplication of matrices `M1^(t)` and `Y`, i.e., `M_1^(t)Y`.
#' @param XtY Input the multiplication of matrices `X^(t)` and `Y`, i.e., `X^(t)Y`.
#' @param penaltyParameterList
#' \itemize{
#'   \item Penalty='PasswayLasso' need two parameters.
#'   \itemize{
#'     \item kappa The L1-norm penalty for pathway Lasso.
#'     \item nu The L2-norm penalty for pathway Lasso.
#'   }
#'   \item Penalty='Network' need one parameter.
#'   \itemize{
#'     \item L The network.
#'   }
#'   \item Penalty='Lasso' don't need other parameters.
#'   \item Penalty='ElasticNet' don't need other parameters.
#' }
#' @return A object, SingleModalityAdmm, with three elements.
#' \itemize{
#'   \item gamma: estimated direct effect.
#'   \item alpha: estimate effect between mediator and independent variables.
#'   \item beta : estimate effect between mediator and dependent variables.
#' }
#' @param verbose A logical value to specify whether to print the iteration process.
#' @param debug A logical value to specify whether to print more iteration processes.
#' @references
#' \enumerate{
#'   \item Zhao, Y., & Luo, X. (2022). Pathway Lasso: pathway estimation and selection with high-dimensional mediators. Statistics and its interface, 15(1), 39.
#' }
#' @examples
#' ## Generate Empirical Data
#' simuData <- One_Modality_Mediation_Data_Gen_fixNonZero(
#'   n = 50, p = 100,
#'   parameters = list(
#'     sigmaY = 1,
#'     sigmaM1 = diag(p),
#'     sizeNonZero = c(3, 3, 4),
#'     alphaMean = c(6, 4, 2),
#'     alphaSd = c(0.1, 0.1, 0.1),
#'     tauMean = c(6,4,2),
#'     tauSd = c(0.1, 0.1, 0.1)
#'   ),
#'   trueValue = list(gamma = 3),
#'   seed = 20231201
#' )
#'
#'
#' ## Parameter Estimation
#' model <- singleModalityAdmm(
#'   X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
#'   rho = 1, lambda1g = 2, lambda1a = 1, lambda1b = 0.1, lambda2a = 1, lambda2b = 1,
#'   penalty = "PathywayLasso", penaltyParameterList = list(kappa = 1, nu = 2),
#'   SIS = FALSE, SIS_thre = 2
#' )
#' @export
singleModalityAdmm <- function(
    X, Y, M1,
    rho=1, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b,
    penalty = "Network", penaltyParameterList,
    maxIter=3000, tol=1e-4,
    verbose = FALSE, debug = FALSE, useCpp = FALSE
) {
  YY <- scale(Y)
  Y.center <- attr(YY,"scaled:center")
  Y.scale <- attr(YY,"scaled:scale")
  XX <- scale(X)
  X.center <- attr(XX,"scaled:center")
  X.scale <- attr(XX,"scaled:scale")

  MM1 <- scale(M1)
  M1.center <- attr(MM1,"scaled:center")
  M1.scale <- attr(MM1,"scaled:scale")

  n <- length(M1)
  p <- ncol(M1)

  preCalcValues <- list(
    XtX = fMatTransProd(XX, XX),
    XtM1 = fMatTransProd(XX, MM1),
    M1tM1PlusRhoInv = fMatInv(fMatTransProd(MM1, MM1) + rho*diag(p)),
    M1tY = fMatTransProd(MM1, YY),
    XtY = fMatTransProd(XX, YY)
  )

  gammaEst <- matrix(coef(lm(YY~XX+MM1))[2:(2+ncol(X)-1)], ncol=1)
  alphaEst <- coef(lm(MM1~XX))[2, , drop=FALSE]
  betaEst <- matrix(sapply(1:p, function(i) coef(lm(YY~XX+MM1[,i]))[3]), ncol=1)
  tauAlphaEst <- matrix(rep(0, p), nrow=1)
  tauBetaEst <- matrix(rep(0, p), ncol=1)


  if (penalty == "Network") {
    if (!("L" %in% names(penaltyParameterList))) {
      stop("penaltyParameterList should contains L for Network penalty")
    }
  } else if (penalty == "PasswayLasso") {
    if (!("kappa" %in% names(penaltyParameterList)) || !("nu" %in% names(penaltyParameterList))) {
      stop("penaltyParameterList should contains kappa and nu for PasswayLasso penalty")
    }
  } else if (penalty == "Lasso") {
    # do nothing
  } else if (penalty == "ElasticNet") {
    # do nothing
  } else {
    stop("No such penalty.")
  }

  iter <- 0L
  estFunctionName <- sprintf("esstimate%s", penalty)
  commonArgs <- list(
    X=XX, Y=YY, M1=MM1, alpha=alphaEst, beta=betaEst, gamma=gammaEst, tauAlpha=tauAlphaEst, tauBeta=tauBetaEst,
    rho = rho, lambda1a=lambda1a, lambda1b=lambda1b, lambda1g=lambda1g, lambda2a=lambda2a, lambda2b=lambda2b,
    XtX=preCalcValues$XtX, XtM1=preCalcValues$XtM1, M1tM1PlusRhoInv=preCalcValues$M1tM1PlusRhoInv, M1tY=preCalcValues$M1tY, XtY=preCalcValues$XtY
  )
  estRes <- do.call(estFunctionName, c(commonArgs, penaltyParameterList))

  if (verbose || debug) {
    cat(sprintf(
      "Initial values - gamma: %.6f, alpha[1]: %.6f, beta[1]: %.6f, tauAlpha[1]: %.6f, tauBeta[1]: %.6f\n",
      gammaEst, alphaEst[1], betaEst[1], tauAlphaEst[1], tauBetaEst[1]
    ))
    cat(sprintf(
      "Iteration %i - gamma: %.6f, alpha[1]: %.6f, beta[1]: %.6f, tauAlpha[1]: %.6f, tauBeta[1]: %.6f\n",
      iter, estRes$gamma, estRes$alpha[1], estRes$beta[1], estRes$tauAlpha[1], estRes$tauBeta[1]
    ))
  }

  cont <- TRUE
  estResOld <- list(
    alphaStep1 = alphaEst+9999, betaStep2 = betaEst+9999,
    alpha = alphaEst+9999, beta = betaEst+9999, gamma = gammaEst+9999,
    tauAlpha = tauAlphaEst+9999, tauBeta = tauBetaEst+9999
  )

  while (cont) {
    iter <- iter + 1L

    estResOld <- estRes
    commonArgs <- list(
      X=XX, Y=YY, M1=MM1, alpha=estRes$alpha, beta=estRes$beta, gamma=estRes$gamma, tauAlpha=estRes$tauAlpha, tauBeta=estRes$tauBeta,
      rho = rho, lambda1a=lambda1a, lambda1b=lambda1b, lambda1g=lambda1g, lambda2a=lambda2a, lambda2b=lambda2b,
      XtX=preCalcValues$XtX, XtM1=preCalcValues$XtM1, M1tM1PlusRhoInv=preCalcValues$M1tM1PlusRhoInv, M1tY=preCalcValues$M1tY, XtY=preCalcValues$XtY
    )
    estRes <- do.call(estFunctionName, c(commonArgs, penaltyParameterList))

    convCondNums <- c(
      sum((estRes$gamma - estResOld$gamma)^2),
      sum((estRes$alphaStep1 - estRes$alpha)^2),
      sum((estRes$alphaStep1 - estResOld$alphaStep1)^2),
      sum((estRes$betaStep2 - estRes$beta)^2),
      sum((estRes$betaStep2 - estResOld$betaStep2)^2)
    )
    cont <- (iter < maxIter) && any(convCondNums >= tol)

    if (verbose || debug) {
      if ((verbose && (iter %% 10 == 0)) || debug) {
        cat(sprintf(
          "Iteration %i - isConv: %i, gamma: %.6f, alpha[1]: %.6f, beta[1]: %.6f, alphaStep1[1]: %.6f, betaStep2[1]: %.6f, tauAlpha[1]: %.6f, tauBeta[1]: %.6f\n",
          iter, !cont, estRes$gamma, estRes$alpha[1], estRes$beta[1], estRes$alphaStep1[1], estRes$betaStep2[1], estRes$tauAlpha[1], estRes$tauBeta[1]
        ))

        if (debug) {
          cat(sprintf(
            "         - converging condition SSE: %.6f, %.6f, %.6f, %.6f, %.6f\n",
            convCondNums[1], convCondNums[2], convCondNums[3], convCondNums[4], convCondNums[5]
          ))
        }
      }
    }
  }

  if (iter >= maxIter) {
    warning("Method does not converge!")
  }

  out <- list(
    alpha = estRes$alpha * M1.scale / X.scale,
    beta = estRes$beta * Y.scale/M1.scale,
    gamma = estRes$gamma * Y.scale/X.scale,
    tauAlpha = estRes$tauAlpha,
    tauBeta = estRes$tauBeta,
    isConv = iter < maxIter,
    niter = iter
  )

  out$interceptAlpha <- mean(M1) - fMatProd(matrix(X.center, nrow=1), out$alpha) * M1.scale / X.scale
  out$interceptBeta <- Y.center - as.numeric(fMatAdd(fMatProd(matrix(X.center, nrow=1), out$gamma), fMatProd(fMatProd(matrix(X.center, nrow=1), out$alpha), out$beta)))
  out$loglik <- getLogLikelihood(X, Y, M1, out$alpha, out$beta, matrix(out$gamma, 1, 1), out$interceptAlpha, out$interceptBeta)
  class(out) <- "SingleModalityAdmm"
  return(out)
}
