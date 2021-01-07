test_that("The number of cores for RcppParallel", {
  expect_error(set_cores("s"),
               "Please enter valid type - but got character")
  expect_error(set_cores(-1),
               "The number of cores is not greater than 1 - but got -1")
  expect_true(set_cores(3))
})
