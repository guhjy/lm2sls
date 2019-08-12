m <- lm2sls(Q ~ P + D, instruments = ~ D + F + A, data=Kmenta)
m1 <- update(m, subset = -1)
test_that("dfbeta computed correctly 1", {
  expect_equal(dfbeta(m)[1, ], coef(m) - coef(m1))
})

m.2sls <- lm2sls(Q ~ P + D, instruments = ~ P + D, data=Kmenta) # OLS
m.lm <- lm(Q ~ P + D, data=Kmenta)

test_that("hatvalues computed correctly", {
  expect_equal(as.vector(hatvalues(m.2sls)), as.vector(hatvalues(m.lm)))
})

test_that("Cook's distances computed correctly", {
  expect_equal(cooks.distance(m.2sls), cooks.distance(m.lm))
})

test_that("dfbeta computed correctly 2", {
  expect_equal(dfbeta(m.2sls), dfbeta(m.lm))
})

m.2sls.w <- lm2sls(Q ~ P + D, instruments = ~ P + D, weights=Q, data=Kmenta) # WLS
m.lm.w <- lm(Q ~ P + D, data=Kmenta, weights=Q)

test_that("hatvalues computed correctly with weights", {
  expect_equal(as.vector(hatvalues(m.2sls.w)), as.vector(hatvalues(m.lm.w)))
})

test_that("Cook's distances computed correctly with weights", {
  expect_equal(cooks.distance(m.2sls.w), cooks.distance(m.lm.w))
})

test_that("rstudent computed correctly", {
  expect_equal(rstudent(m.2sls), rstudent(m.lm))
})

mw <- lm2sls(Q ~ P + D, instruments = ~ D + F + A, weights=Q, data=Kmenta)
m1w <- update(mw, subset = -10)
test_that("dfbeta computed correctly with weights", {
  expect_equal(dfbeta(mw)[10, ], coef(mw) - coef(m1w))
})

test_that("rstudent computed correctly with weights", {
  expect_equal(rstudent(m.2sls.w), rstudent(m.lm.w))
})

