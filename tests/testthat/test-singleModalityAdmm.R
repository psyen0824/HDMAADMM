context("test - singleModalityAdmm")

## Generate Empirical Data
simuData <- modalityMediationDataGen(
  n = 10, p = 5,
  sizeNonZero = c(1),
  alphaMean = c(6),
  betaMean = c(6),
  seed = 20231201
)

modelElasticNet <- singleModalityAdmm(
  X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
  rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
  penalty = "ElasticNet"
)

test_that("singleModalityAdmm - EN", {
  expect_equal(modelElasticNet$alpha, matrix(c(5.56137712, 4.03423839, 1.29329287, 0.16793167, 0.0), nrow=1))
  expect_equal(modelElasticNet$beta, matrix(c(4.13741388, 1.45642770, 7.57421846, 0.0, 0.0), ncol=1))
  expect_equal(modelElasticNet$gamma, matrix(0, nrow=1, ncol=1))
  expect_equal(predict(modelElasticNet, matrix(c(0, 1))), matrix(c(1.90794012, 40.58891832), ncol=1))
  expect_equal(class(modelElasticNet), "SingleModalityAdmm")
})

test_that("singleModalityAdmm - Expect Error", {
  expect_error(
    singleModalityAdmm(
      X = matrix(simuData$MediData$X, ncol=2,nrow=10), Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
      rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      penalty = "ElasticNet"
    )
  )

  expect_error(
    singleModalityAdmm(
      X = simuData$MediData$X[1:9, ], Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
      rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      penalty = "ElasticNet"
    )
  )

  expect_error(
    singleModalityAdmm(
      X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1[1:9, ],
      rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      penalty = "ElasticNet"
    )
  )

  expect_error(
    singleModalityAdmm(
      X = simuData$MediData$X, Y = simuData$MediData$Y[1:9, ], M1 = simuData$MediData$M1,
      rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      penalty = "ElasticNet"
    )
  )

  expect_error(
    singleModalityAdmm(
      X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
      rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      penalty = "ElasticNet1"
    )
  )

  expect_error(
    singleModalityAdmm(
      X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
      rho = 1, lambda1a = 1, lambda1b = NA_real_, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      penalty = "ElasticNet"
    )
  )

  expect_error(
    singleModalityAdmm(
      X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
      rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      maxIter = NA_integer_, penalty = "ElasticNet"
    )
  )

  expect_error(
    singleModalityAdmm(
      X = simuData$MediData$X, Y = simuData$MediData$Y, M1 = simuData$MediData$M1,
      rho = 1, lambda1a = 1, lambda1b = 0.5, lambda1g = 2, lambda2a = 1, lambda2b = 1,
      penalty = "ElasticNet",
      verboseOptions = list(numIter = NA_integer_)
    )
  )
})
