# generate 1-D data with a given seed
set.seed(1234)
originalPar <- par(no.readonly = TRUE)
num_cores <- 2

p <- q <- 20
n <- 100
x1 <- matrix(seq(-7, 7, length = p), nrow = p, ncol = 1)
x2 <- matrix(seq(-7, 7, length = q), nrow = q, ncol = 1)
u <- exp(-x1^2) / norm(exp(-x1^2), "F")
v <- exp(-(x2 - 2)^2) / norm(exp(-(x2 - 2)^2), "F")
Sigma <- array(0, c(p + q, p + q))
Sigma[1:p, 1:p] <- diag(p)
Sigma[(p + 1):(p + q), (p + 1):(p + q)] <- diag(p)
Sigma[1:p, (p + 1):(p + q)] <- u %*% t(v)
Sigma[(p + 1):(p + q), 1:p] <- t(Sigma[1:p, (p + 1):(p + q)])
noise <-
  MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = 0.001 * diag(p + q))
Y <- MASS:::mvrnorm(n, mu = rep(0, p + q), Sigma = Sigma) + noise
Y1 <- Y[, 1:p]
Y2 <- Y[, -(1:p)]

cv_1D <- spatmca(
  x1,
  x2,
  Y1,
  Y2,
  K = 1,
  num_cores = num_cores
)
usedNumberCores <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))

# Test the result
tol <- 1e-4
test_that("Selected tuning parameters", {
  expect_lte(abs(1 - norm(cv_1D$Uestfn, "F")), tol)
  expect_lte(abs(1 - norm(cv_1D$Vestfn, "F")), tol)
  expect_lte(abs(cv_1D$Khat - 1), tol)
  expect_null(cv_1D$cvall)
})

test_that("Number of threads", {
  expect_equal(num_cores, usedNumberCores)
})


test_that("CV plot", {
  expect_error(
    plot.spatmca("test"),
    cat("Invalid object! Please enter a `spatmca` object")
  )
  expect_equal(class(plot.spatmca(cv_1D)), "list")
})


test_that("spatmca handles input with different dimensions", {
  x1_diff_dim <- matrix(1:10, nrow = 2, ncol = 5)  # Adjust dimensions to be different
  expect_error(
    spatmca(x1_diff_dim, x2, Y1, Y2),
    "The number of rows of x1 should be equal to the number of columns of Y1."
  )
})


test_that("Setting tau2u and tau2v when both are provided", {
  tau2u <- c(0.1, 0.5, 1.0)
  tau2v <- c(0.1, 0.5, 1.0)
  result <- spatmca(x1, x2, Y1, Y2, tau2u = tau2u, tau2v = tau2v)
  ntau1u <- length(result$tau1u)
  ntau1v <- length(result$tau1v)
  ntau2u <- length(result$tau2u)
  ntau2v <- length(result$tau2v)
  expect_equal(ntau1u, 11)
  expect_equal(ntau1v, 11)  
  expect_equal(ntau2u, 3)
  expect_equal(ntau2v, 3)
})
test_that("Setting tau1u and tau1v when both are provided", {
  tau1u <- c(0.1, 0.5, 1.0)
  tau1v <- c(0.1, 0.5, 1.0)
  result <- spatmca(x1, x2, Y1, Y2, tau1u = tau1u, tau1v = tau1v)
  ntau1u <- length(result$tau1u)
  ntau1v <- length(result$tau1v)
  ntau2u <- length(result$tau2u)
  ntau2v <- length(result$tau2v)
  expect_equal(ntau1u, 3)
  expect_equal(ntau1v, 3)  
  expect_equal(ntau2u, 11)
  expect_equal(ntau2v, 11)
})

test_that("Setting tau1u when it is NULL", {
  result <- spatmca(x1, x2, Y1, Y2, tau1u = c(0.1, 0.5, 1.0), M = 3)
  expect_equal(length(result$tau1u), 3)
})

test_that("Setting tau1v when it is NULL", {
  result <- spatmca(x1, x2, Y1, Y2, tau1v = c(0.1, 0.5, 1.0), M = 3)
  expect_equal(length(result$tau1v), 3)
})

test_that("Setting tau2u when it is NULL", {
  result <- spatmca(x1, x2, Y1, Y2, tau2u = c(0.1, 0.5, 1.0), M = 3)
  expect_equal(length(result$tau2u), 3)
})

test_that("Setting tau2v when it is NULL", {
  result <- spatmca(x1, x2, Y1, Y2, tau2v = c(0.1, 0.5, 1.0), M = 3)
  expect_equal(length(result$tau2v), 3)
})
