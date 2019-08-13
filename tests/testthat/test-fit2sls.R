
m.2sls <- with(Kmenta, fit2sls(Q, cbind(1, P, D), cbind(1, D, F, A)))
m.stage1 <- lm(cbind(1, P, D) ~ cbind(1, D, F, A) - 1, data=Kmenta)
m.stage2 <- lm(Kmenta$Q ~ fitted(m.stage1) - 1)


test_that("2SLS coefficients are computed correctly", {
  expect_equal(as.vector(coef(m.2sls)), as.vector(coef(m.stage2)))
})

se.2sls <- as.vector(sqrt(diag(m.2sls$vcov)))
se <- as.vector(sqrt(diag(vcov(m.stage2))))
residuals <- as.vector(with(Kmenta, Q - cbind(1, P, D) %*% coef(m.stage2)))
se <- sqrt(sum(residuals^2)/(nrow(Kmenta) - 3))*se/sigma(m.stage2)

test_that("2SLS coefficient standard errors are computed correctly", {
  expect_equal(se.2sls, se)
})

test_that("model residuals are computed correctly", {
  expect_equal(residuals, residuals(m.2sls))
})


m.2sls.w <- with(Kmenta, fit2sls(Q, cbind(1, P, D), cbind(1, D, F, A), wt=Q))
m.stage1.w <- lm(cbind(1, P, D) ~ cbind(1, D, F, A) - 1, weights=Q, data=Kmenta)
m.stage2.w <- lm(Kmenta$Q ~ fitted(m.stage1.w) - 1, weights=Kmenta$Q)

test_that("2SLS coefficients are computed correctly with weights", {
  expect_equal(as.vector(coef(m.2sls.w)), as.vector(coef(m.stage2.w)))
})

se.2sls.w <- as.vector(sqrt(diag(m.2sls.w$vcov)))
se.w <- as.vector(sqrt(diag(vcov(m.stage2.w))))
residuals <- as.vector(with(Kmenta, Q - cbind(1, P, D) %*% coef(m.stage2.w)))
se.w <- sqrt(sum(Kmenta$Q*residuals^2)/(nrow(Kmenta) - 3))*se.w/sigma(m.stage2.w)

test_that("2SLS coefficient standard errors are computed correctly with weights", {
  expect_equal(se.2sls.w, se.w)
})
