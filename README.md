HDMAADMM
====

Parameter Estimation of High-dimensional Mediation Models Based on Alternating Direction Method of Multipliers (ADMM)

Installation
------------

You can install:

-   the latest development version from github with

``` r
install.packages("remotes")
remotes::install_github("psyen0824/HDMAADMM")
```

If you encounter a bug, please file a reproducible example on [github](https://github.com/psyen0824/HDMAADMM/issues).

#### If you would like to have `oneMKL` support. Please use below commands to install:

``` r
install.packages(
  c("oneMKL", "oneMKL.MatrixCal"), 
  repos=c("https://cloud.r-project.org/", "https://R-OneMKL.github.io/drat"), 
  type = "source"
)
remotes::install_github("psyen0824/HDMAADMM", ref = "oneMKL")
```

examples
--------

``` r
simuData <- modalityMediationDataGen(seed = 20231201)

modelElasticNet <- singleModalityAdmm(
  X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
  rho = 1, lambda1a = 1, lambda1b = 0.1, lambda1g = 2, lambda2a = 1, lambda2b = 1,
  penalty = "ElasticNet"
)

fitted(modelElasticNet)                        # fitted values
predict(modelElasticNet, matrix(c(0, 1), ncol = 1))  # predict values
```
