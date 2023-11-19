library(oneMKL)
library(oneMKL.MatrixCal)
library(Rcpp)

softThreshold <- function(b, lambda) sign(b) * ifelse(abs(b) <= lambda, 0, abs(b)-lambda)

getLogLikelihood <- function(X, Y, M1, alpha, beta, gamma, interceptAlpha, interceptBeta) {
  # A <- fMatSubtract(M1, fMatProd(X, alpha))
  A <- fMatSubtract(M1, sweep(fMatProd(X, alpha), 2, interceptAlpha, `+`))
  l1 <- 0
  for (i in seq_len(ncol(A))) {
    l1 <- l1 + (-1/2) * crossprod(A[ , i])
  }

  # B <-  fMatSubtract(Y, fMatAdd(fMatProd(X, gamma), fMatProd(M1, beta)))
  B <-  fMatSubtract(Y, fMatAdd(fMatProd(X, gamma), fMatProd(M1, beta))) - interceptBeta[1]
  l2 <- (-1/2)*fMatTransProd(B, B)

  return(list(l=as.numeric(l1+l2), l1=as.numeric(l1), l2=as.numeric(l2)))
}

rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean))) {
  R <- chol(sigma, pivot = TRUE)
  R[, order(attr(R, "pivot"))]

  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% R
  retval <- sweep(retval, 2, mean, "+")
  retval
}
