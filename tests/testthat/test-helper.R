defaultNumber <- RcppParallel::defaultNumThreads()
tol <- 1e-6

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

# Test invalid input
x1 <- matrix(seq(-7, 7, length = 4), nrow = 4, ncol = 1)
x2 <- matrix(seq(-7, 7, length = 5), nrow = 5, ncol = 1)
Y1 <- matrix(rnorm(n = 100 * 4), 100, 4)
Y2 <- matrix(rnorm(n = 100 * 5), 100, 5)
M <- 5
test_that("check input of spatmca", {
  expect_null(checkInputData(x1, x2, Y1, Y2, M))
  expect_error(
    checkInputData(matrix(1:10, nrow = 1), x2, Y1, Y2, M),
    cat("The number of rows of x1 should be equal to the number of columns of Y1.")
  )
  expect_error(
    checkInputData(x1, matrix(1:10, nrow = 1), Y1, Y2, M),
    cat("The number of rows of x2 should be equal to the number of columns of Y2.")
  )
  expect_error(
    checkInputData(matrix(1:10, nrow = 2), x2, matrix(rnorm(n = 100 * 2), 100, 2), Y2, M),
    cat("Number of locations must be larger than 2.")
  )
  expect_error(
    checkInputData(matrix(1:20, nrow = 4), x2, matrix(rnorm(n = 100 * 4), 100, 4), Y2, M),
    cat("Dimension of locations must be less 4.")
  )
  expect_error(
    checkInputData(x1, x2, matrix(rnorm(n = 50 * 4), 50, 4), Y2, M),
    cat("The numbers of sample sizes of both data should be equal")
  )
  expect_error(
    checkInputData(x1, x2, Y1, Y2, 101),
    cat("Number of folds must be less than sample size, but got M = 101")
  )
})

# Test detrend
test_that("check detrending", {
  expect_equal(detrend(Y1, FALSE), Y1)
  expect_lte(sum(detrend(Y1, TRUE)), tol)
})

# Test plot functions
set.seed(1234)
originalPar <- par(no.readonly = TRUE)
num_cores <- 2

p <- q <- 4
n <- 20
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
cv_data <- expand.grid(u = cv_1D$tau1u, v = cv_1D$tau1v)
cv_data$cv <- as.vector(cv_1D$cv1)

plot_obj <- plot_cv_field(cv_data, "first variate")
test_that("GGplot objects", {
  expect_true("ggplot" %in% class(plot_obj))
})

plot_sequentially(list(plot_obj))
newPar <- par(no.readonly = TRUE)
test_that("Envirorment setting", {
  expect_equal(originalPar, newPar)
})

