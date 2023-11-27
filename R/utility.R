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

#' Helper function to convert Weight Matrix to Laplacian Matrix
#'
#' @param W The weight matrix for n nodes which should be \code{n}x\code{n} matrix.
#' @return L \code{n}x\code{n} Laplacian matrix.

WeightToLaplacian <- function(W) {
  if ((nrow(W) != ncol(W)) || !isSymmetric(W)) {
    stop("W should be a square symmetric matrix.")
  }
  d <- colSums(W)
  L <- -W / sqrt(d %o% d)
  diag(L) <- diag(L) + 1.0
  return(L)
}

getLogLikelihood <- function(X, Y, M1, alpha, beta, gamma) {
  A <- M1 - fMatProd(X, alpha)
  l1 <-  (-1/2) * sum(diag(fMatTransProd(A, A)))
  B <- Y - fMatProd(X, gamma) - fMatProd(M1, beta)
  l2 <- (-1/2)*fMatTransProd(B, B)
  return(list(l=as.numeric(l1+l2), l1=as.numeric(l1), l2=as.numeric(l2)))
}
