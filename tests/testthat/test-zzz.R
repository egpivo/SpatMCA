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
