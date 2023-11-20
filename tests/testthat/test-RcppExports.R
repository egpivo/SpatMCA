library(testthat)
library(Rcpp)
library(RcppArmadillo)

# Create a unit test
test_that("spatmcacvall_rcpp returns the correct result", {
  # Set up test data (replace this with your actual test data)
  sxr <- matrix(1:20, nrow = 10, ncol = 2)
  syr <- matrix(21:40, nrow = 10, ncol = 2)
  Xr <- matrix(rnorm(100), nrow = 10, ncol = 10)
  Yr <- matrix(rnorm(100), nrow = 10, ncol = 10)
  M <- 5
  K <- 2
  tau1ur <- c(0.1, 0.5, 1.0)
  tau2ur <- c(0.2, 0.8)
  tau1vr <- c(0.1, 0.5, 1.0)
  tau2vr <- c(0.2, 0.8)
  nkr <- c(5, 10, 15)
  maxit <- 100
  tol <- 1e-6
  l2ur <- c(0.1, 0.5, 1.0)
  l2vr <- c(0.2, 0.8)
  
  # Call the Rcpp function
  result <- spatmcacvall_rcpp(sxr, syr, Xr, Yr, M, K, tau1ur, tau2ur, tau1vr, tau2vr, nkr, maxit, tol, l2ur, l2vr)
  
  # Check if the result is a list with expected components
  expect_is(result, "list")
  expect_true("cvall" %in% names(result))
  expect_true("Uest" %in% names(result))
  expect_true("Vest" %in% names(result))
  expect_true("Dest" %in% names(result))
  expect_true("cvtau1u" %in% names(result))
  expect_true("cvtau2u" %in% names(result))
  expect_true("cvtau1v" %in% names(result))
  expect_true("cvtau2v" %in% names(result))
  
  # Add more specific checks based on the expected output
  # ...
  
})
