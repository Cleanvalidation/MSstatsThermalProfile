test_that("Correlation to Covariance Matrix works", {
  R=matrix(c(1,0.5,0.8,0.5,1,0.3,0.8,0.3,1),nrow=3,ncol=3,byrow=TRUE)
  v = c(4,2,1)

  S_expected = matrix(c(4,sqrt(2),1.6,sqrt(2),2,0.3*sqrt(2),1.6,0.3*sqrt(2),1),nrow=3,ncol=3,byrow=TRUE)
  S_actual = cor2cov(R,v)
  expect_equal(S_actual, S_expected)
})
