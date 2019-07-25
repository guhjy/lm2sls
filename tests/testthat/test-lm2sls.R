m.fit2sls <- with(Kmenta, fit2sls(Q, cbind(1, P, D), cbind(1, D, F, A)))
m.lm2sls <- lm2sls(Q ~ P + D, instruments = ~ D + F + A, data=Kmenta)

test_that("2SLS coefficients are computed correctly", {
  expect_equal(as.vector(coef(m.fit2sls)), as.vector(coef(m.lm2sls)))
})

test_that("2SLS coefficient covariance matrix is computed correctly", {
  expect_equal(as.vector(m.fit2sls$vcov), as.vector(vcov(m.lm2sls)))
})
