defaultNumber <- RcppParallel::defaultNumThreads()
test_that("The number of cores for RcppParallel", {
  expect_error(
    set_cores("s"),
    "Please enter valid type - but got character"
  )
  expect_error(
    set_cores(-1),
    "The number of cores is not greater than 1 - but got -1"
  )
  expect_error(
    set_cores(defaultNumber + 1),
    cat("The input number of cores is invalid - default is ", defaultNumber)
  )
  expect_true(set_cores(defaultNumber))
  expect_null(set_cores())
})
