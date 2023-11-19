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
##
## You should have received a copy of the GNU General Public License
## along with HDMAADMM If not, see <http://www.gnu.org/licenses/>.

#' `HDMAADMM` Package
#'
#' This package enables the estimation of single-modality high-dimensional mediation models.
#' We employ penalized maximum likelihood and solve the estimation using the
#' Alternating Direction Method of Multipliers (ADMM) to provide high-dimensional mediator estimates.
#' To improve the sensitivity and specificity of non-zero mediators,
#' we offer the sure independence screening (SIS) function for dimension reduction.
#' The available penalty options include Lasso, Elastic Net, Pathway Lasso, and Network-constrained Penalty.
#' The methods employed in the package are based on Lasso (doi.org/10.1111/j.2517-6161.1996.tb02080.x),
#' Elastic Net (doi.org/10.1111/j.1467-9868.2005.00503.x), Pathway Lasso (10.4310/21-sii673),
#' Network Penalty (doi.org/10.1093/bioinformatics/btn081), ADMM (doi.org/10.1561/2200000016),
#' and SIS (doi.org/10.1111/j.1467-9868.2008.00674.x).
#'
#' @docType package
#' @name HDMAADMM-package
#' @useDynLib oneMKL.MatrixCal
#' @importFrom Rcpp evalCpp
#' @importFrom RcppEigen RcppEigen.package.skeleton
NULL
