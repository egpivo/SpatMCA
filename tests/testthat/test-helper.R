defaultNumber <- RcppParallel::defaultNumThreads()
test_that("The number of cores for RcppParallel", {
  expect_error(
    setCores("s"),
    "Please enter valid type - but got character"
  )
  expect_error(
    setCores(-1),
    "The number of cores is not greater than 1 - but got -1"
  )
  expect_error(
    setCores(defaultNumber + 1),
    cat("The input number of cores is invalid - default is ", defaultNumber)
  )
  expect_true(setCores(defaultNumber))
  expect_null(setCores())
})
