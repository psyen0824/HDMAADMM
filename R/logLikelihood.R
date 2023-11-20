getLogLikelihood <- function(X, Y, M1, alpha, beta, gamma, interceptAlpha, interceptBeta) {
  A <- M1 - sweep(fMatProd(X, alpha), 2, interceptAlpha, `+`)
  l1 <- 0
  for (i in seq_len(ncol(A))) {
    l1 <- l1 + (-1/2) * crossprod(A[ , i])
  }
  B <-  Y - fMatProd(X, gamma) + fMatProd(M1, beta) - interceptBeta[1]
  l2 <- (-1/2)*fMatTransProd(B, B)

  return(list(l=as.numeric(l1+l2), l1=as.numeric(l1), l2=as.numeric(l2)))
}
