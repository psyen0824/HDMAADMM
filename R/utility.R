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

#' Helper function to convert Adjacency Matrix to Laplacian Matrix
#'
#' @param A The Adjacency matrix for n nodes which should be \code{n}x\code{n} matrix.
#' @return A \code{n}x\code{n} Laplacian matrix.
#' @examples
#' (A <- matrix(c(1, 1.5, 0, 1.5, 2, 0, 0, 0, 3), 3, 3))
#' (L <- adjacencyToLaplacian(A))
adjacencyToLaplacian <- function(A) {
  if ((nrow(A) != ncol(A)) || !isSymmetric(A)) {
    stop("A should be a square symmetric matrix.")
  }
  d <- colSums(A)
  L <- -A / sqrt(d %o% d)
  diag(L) <- 1.0
  return(L)
}
