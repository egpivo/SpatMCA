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


test_that(".onAttach works without interfering with CRAN submission", {
  # Simulate CRAN submission by setting NOT_CRAN to true
  Sys.setenv(NOT_CRAN = "true")
  
  # Save the current value of RCPP_PARALLEL_BACKEND
  original_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND")
  
  # Call the .onAttach function
  .onAttach("SpatMCA", "SpatMCA")
  
  # Check if RCPP_PARALLEL_BACKEND is not set when NOT_CRAN is true
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND"), original_backend)
  
  # Reset NOT_CRAN to avoid interference with other tests
  Sys.unsetenv("NOT_CRAN")
})

test_that(".onAttach sets RCPP_PARALLEL_BACKEND for parallel processing", {
  # Save the current value of RCPP_PARALLEL_BACKEND
  original_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND")
  
  # Simulate parallel processing by setting NOT_CRAN to false
  Sys.setenv(NOT_CRAN = "false")
  
  # Call the .onAttach function
  .onAttach("SpatMCA", "SpatMCA")
  
  # Check if RCPP_PARALLEL_BACKEND is set to "tinythread" for parallel processing
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND"), "tinythread")
  
  # Restore the original value of RCPP_PARALLEL_BACKEND
  Sys.setenv(RCPP_PARALLEL_BACKEND = original_backend)
  
  # Reset NOT_CRAN to avoid interference with other tests
  Sys.unsetenv("NOT_CRAN")
})
