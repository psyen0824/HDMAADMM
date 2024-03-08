context("test - utility")

weightMat <- matrix(c(1,2,3,2,4,5,3,5,6), nrow = 3, ncol = 3)
expectedLaplacianMat <- matrix(c(0.83333333, -0.24618298, -0.32732684, -0.24618298, 0.63636364, -0.40291148, -0.32732684, -0.40291148, 0.57142857), nrow = 3, ncol = 3)
test_that("weightToLaplacian", {
  expect_equal(weightToLaplacian(weightMat), expectedLaplacianMat)
  expect_error(weightToLaplacian(matrix(1:9, 3)))
  expect_error(weightToLaplacian(matrix(1:12, 4)))
})

test_that("checkPenaltyParameterList", {
  expect_error(HDMAADMM:::checkPenaltyParameterList(list(x = 4), "y", "test"))
  expect_error(HDMAADMM:::checkPenaltyParameterList(list(x = NA), "x", "test"))
  expect_error(HDMAADMM:::checkPenaltyParameterList(list(x = Inf), "x", "test"))
})
