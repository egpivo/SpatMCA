library(testthat)
library(SpatMCA)

test_that("onAttach sets RCPP_PARALLEL_BACKEND to tinythread on CRAN", {
  old_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND")
  Sys.setenv(NOT_CRAN = "true")
  Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
  
  # Load the package, triggering .onAttach
  library(SpatMCA)
  
  # Check that RCPP_PARALLEL_BACKEND is set to tinythread when NOT_CRAN is "true"
  if (Sys.getenv("NOT_CRAN") == "true") {
    expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND"), "tinythread")
  } else {
    # If NOT_CRAN is not "true", RCPP_PARALLEL_BACKEND should not be set
    expect_true(is.null(Sys.getenv("RCPP_PARALLEL_BACKEND")))
  }
  
  # Reset the environment variable to its original value
  Sys.setenv(NOT_CRAN = "false")
  Sys.setenv(RCPP_PARALLEL_BACKEND = old_backend)
})

test_that("Setting RCPP_PARALLEL_BACKEND in non-CRAN environment", {
  # Mock NOT_CRAN environment variable
  mock_not_cran <- "false"
  old_env <- Sys.getenv("NOT_CRAN")
  on.exit(Sys.setenv(NOT_CRAN = old_env), add = TRUE)
  Sys.setenv(NOT_CRAN = mock_not_cran)
  
  # Call the code from zzz.R that sets RCPP_PARALLEL_BACKEND
  source("../../R/zzz.R", local = TRUE)
  
  # Print the current value of RCPP_PARALLEL_BACKEND for debugging
  print(Sys.getenv("RCPP_PARALLEL_BACKEND"))
  
  # Check if RCPP_PARALLEL_BACKEND is set to "tinythread"
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND"), "tinythread")
})


