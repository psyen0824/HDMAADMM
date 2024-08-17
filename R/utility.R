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

checkPenaltyParameterList <- function(ppl, parameterNames, penaltyName) {
  for (pn in parameterNames) {
    if (!(pn %in% names(ppl))) {
      stop(sprintf("penaltyParameterList should contains %s for %s penalty!", pn, penaltyName))
    }
    if (any(is.na(ppl[[pn]]) | is.infinite(ppl[[pn]]))) {
      stop(sprintf("penaltyParameterList$%s should be finite non-nan numeric for %s penalty!", pn, penaltyName))
    }
  }
}

