#' Helper function to convert Weight Matrix to Laplacian Matrix
#'
#' @param W The weight matrix for n nodes which should be \code{n}x\code{n} matrix.
#' @return L \code{n}x\code{n} Laplacian matrix.
#' @examples
#' set.seed(20231201)
#' p <- 5
#' W <- matrix(0, nrow = p, ncol = p)
#' W[lower.tri(W)] <- runif(p*(p-1)/2, 0, 1)
#' W[upper.tri(W)] <- t(W)[upper.tri(W)]
#' diag(W) <- 1
#' (L <- weightToLaplacian(W))
#' @export
weightToLaplacian <- function(W) {
  if ((nrow(W) != ncol(W)) || !isSymmetric(W)) {
    stop("W should be a square symmetric matrix.")
  }
  d <- colSums(W)
  L <- -W / sqrt(d %o% d)
  diag(L) <- diag(L) + 1.0
  return(L)
}

# from fmsb
nagelkerkeRSequared <- function(obj) {
  stopifnot("lm" %in% class(obj))
  n <- nrow(obj$model)
  return((1-exp((obj$dev-obj$null)/n))/(1-exp(-obj$null/n)))
}

#' Function Generate Laplacian Matrix
#'
#' @param X,Y,M1 The input data for the single modality mediation model. Details see \code{\link{singleModalityAdmm}}.
#' @param type A string to specify the generated Laplacian matrix is for \eqn{\alpha} or \eqn{\beta}.
#' @export
#' @importFrom stats lm glm
generateLaplacianMatrix <- function(X, Y, M1, type = "beta") {
  p <- ncol(M1)
  W <- matrix(0, p, p)
  if (type == "beta") {
    for (i in 1:p) {
      mdl <- lm(Y ~ M1[ ,i])
      W[i, i] <- 1 - summary(mdl)$coefficients[2, 4]
    }
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        mdl <- lm(Y ~ M1[ ,i] * M1[ ,j])
        W[i, j] <- 1 - summary(mdl)$coefficients[4, 4]
      }
    }
    W[lower.tri(W)] <- t(W)[lower.tri(W)]
    return(weightToLaplacian(W))
  } else {
    for (i in 1:p) {
      mdl <- glm(X ~ M1[ ,i]  , family = "binomial")
      W[i, i] <- 1 - summary(mdl)$coefficients[2, 4]
    }
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        mdl <- glm(X ~ M1[ ,i] * M1[ ,j], family = "binomial")
        W[i, j] <- 1 - summary(mdl)$coefficients[4, 4]
      }
    }
    W[lower.tri(W)] <- t(W)[lower.tri(W)]
    return(weightToLaplacian(W))
  }
}
