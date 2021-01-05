# generate 1-D data with a given seed
set.seed(1234)

p <- q <- 20
n <- 100
x1 <- matrix(seq(-7, 7, length = p), nrow = p, ncol = 1)
x2 <- matrix(seq(-7, 7, length = q), nrow = q, ncol = 1)
u <- exp(-x1^2)/norm(exp(-x1^2), "F")
v <- exp(-(x2 - 2)^2)/norm(exp(-(x2 - 2)^2), "F")
Sigma <- array(0, c(p + q, p + q))
Sigma[1:p, 1:p] <- diag(p)
Sigma[(p + 1):(p + q), (p + 1):(p + q)] <- diag(p)
Sigma[1:p, (p + 1):(p + q)] <- u%*%t(v)
Sigma[(p + 1):(p + q), 1:p] <- t(Sigma[1:p, (p + 1):(p + q)])
noise <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = 0.001*diag(p + q))
Y <- MASS:::mvrnorm(n, mu = rep(0, p + q), Sigma = Sigma) + noise
Y1 <- Y[,1:p]
Y2 <- Y[,-(1:p)]

cv_1D <- spatmca(x1, x2, Y1, Y2)


# Test the result
tol <- 1e-5
test_that("Selected tuning parameters", {
  expect_lte(abs(cv_1D$Dest - 0.8845911), tol)
  expect_lte(abs(1 - norm(cv_1D$Uestfn, "F")), tol)
  expect_lte(abs(1 - norm(cv_1D$Vestfn, "F")), tol)
  expect_lte(abs(cv_1D$Khat - 1), tol)
  expect_null(cv_1D$cvall)
})
