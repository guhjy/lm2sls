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

test_that("rstudent computed correctly", {
  expect_equal(rstudent(m.2sls), rstudent(m.lm))
})
