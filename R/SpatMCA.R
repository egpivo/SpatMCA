#' @title  Regularized spatial MCA
#'
#' @description Produce spatial coupled patterns at the designated locations according to the specified tuning parameters or the tuning parameters selected by M-fold cross-validation.
#'
#' @param x1 Location matrix (\eqn{p \times d}) correponding to Y1. Each row is a location. \eqn{d=1,2} is the dimension of locations.
#' @param x2 Location matrix (\eqn{q \times d}) correponding to Y2. Each row is a location.
#' @param Y1 Data matrix (\eqn{n \times p}) of the first variable stores the values at \eqn{p} locations with sample size \eqn{n}.
#' @param Y2 Data matrix (\eqn{n \times q}) of the second variable stores the values at \eqn{q} locations with sample size \eqn{n}.
#' @param M Optional number of folds; default is 5.
#' @param K Optional user-supplied number of coupled patterns; default is NULL. If K is NULL or doesSelectRank is TRUE, K is selected automatically.
#' @param doesSelectRank If TRUE, K is selected automatically; otherwise, doesSelectRank is set to be user-supplied K. Default depends on user-supplied K.
#' @param tau1u Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence correponding to Y1. If NULL, 10 tau1u values in a range are used.
#' @param tau2u Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence correponding to Y1. If NULL, 10 tau2u values in a range are used.
#' @param tau1v Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence correponding to Y2. If NULL, 10 tau1v values in a range are used.
#' @param tau2v Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence correponding to Y2. If NULL, 10 tau2v values in a range are used.
#' @param x1new New location matrix correponding to Y1. If NULL, it is x1.
#' @param x2new New location matrix correponding to Y2. If NULL, it is x2.
#' @param center If TRUE, center the columns of Y. Default is FALSE.
#' @param plot.cv If TRUE, plot the cv values. Default is FALSE.
#' @param maxit Maximum number of iterations. Default value is 100.
#' @param thr Threshold for convergence. Default value is \eqn{10^{-4}}.
#' @param doesSelectAllTuningParameters If TRUE, The K-fold CV performs to select 4 tuning parameters simultaneously. Default value is FALSE.
#' @param numCores Number of cores used to parallel computing. Default value is NULL (See `RcppParallel::defaultNumThreads()`)
#'
#' @return A list of objects including 
#' \item{Uestfn}{Estimated patterns for Y1 at the new locations, x1New.}
#' \item{Vestfn}{Estimated patterns for Y2 at the new locations, x2New.}
#' \item{Dest}{Estimated singular values.}
#' \item{crosscov}{Estimated cross-covariance matrix between Y1 and Y2.}
#' \item{stau1u}{Selected tau1u.}
#' \item{stau2u}{Selected tau2u.}
#' \item{stau1v}{Selected tau1v.}
#' \item{stau2v}{Selected tau2v.}
#' \item{cv1}{cv socres for tau1u and tau1v when doesSelectAllTuningParameters is FALSE.}
#' \item{cv2}{cv socres for tau2u and tau2v when doesSelectAllTuningParameters is FALSE.}
#' \item{cvall}{cv socres for tau1u, tau2u, tau1v and tau2v when doesSelectAllTuningParameters is TRUE.}
#' \item{tau1u}{Sequence of tau1u-values used in the process.}
#' \item{tau2u}{Sequence of tau2u-values used in the process.}
#' \item{tau1v}{Sequence of tau1v-values used in the process.}
#' \item{tau2v}{Sequence of tau2v-values used in the process.}
#'
#' @details The optimization problem is
#' \deqn{\max_{\mathbf{U}, \mathbf{V}} \frac{1}{n}\mbox{tr}(\mathbf{U}'\mathbf{Y}'_1\mathbf{Y}_2\mathbf{V}) - \tau_{1u}\mbox{tr}(\mathbf{U}'\mathbf{\Omega}_1\mathbf{U}) - \tau_{2u}\sum_{k=1}^K\sum_{j=1}^{p} |u_{jk}|- \tau_{1v}\mbox{tr}(\mathbf{V}'\mathbf{\Omega}_2\mathbf{V})-\tau_{2v}\sum_{k=1}^K\sum_{j=1}^{q} |v_{jk}|,}
#' \eqn{\mbox{subject to $ \mathbf{U}'\mathbf{U}=\mathbf{V}'\mathbf{V}=\mathbf{I}_K$,}} where \eqn{\mathbf{Y}_1} and \eqn{\mathbf{Y}_2} are two data matrices, \eqn{{\mathbf{\Omega}}_1} and \eqn{{\mathbf{\Omega}}_2} are two smoothness matrix, \eqn{\mathbf{V}=\{v_{jk}\}}, and \eqn{\mathbf{U}=\{u_{jk}\}}.
#' @export
#' @author Wen-Ting Wang and Hsin-Cheng Huang
#' @references Wang, W.-T. and Huang, H.-C. (2017). Regularized principal component analysis for spatial data. \emph{Journal of Computational and Graphical Statistics} \bold{26} 14-25.
#' @examples
#' ###### 1D: regular locations
#' x_1D <- as.matrix(seq(-5, 5, length = 50))
#' Phi_1D <- exp(-x_1D^2) / norm(exp(-x_1D^2), "F")
#' set.seed(1234)
#' Y_1D <- rnorm(n = 100, sd = 3) %*% t(Phi_1D) + matrix(rnorm(n = 100 * 50), 100, 50)
#' cv_1D <- spatpca(x = x_1D, Y = Y_1D, numCores = 2)
#' plot(x_1D, cv_1D$eigenfn[, 1], type = "l", main = "1st eigenfunction")
#' lines(x_1D, svd(Y_1D)$v[, 1], col = "red")
#' legend("topleft", c("SpatPCA", "PCA"), lty = 1:1, col = 1:2)
#' 
#' \dontrun{
#'   ### 1D: artificial irregular locations
#'   rm_loc <- sample(1:50, 20)
#'   x_1Drm <- x_1D[-rm_loc]
#'   Y_1Drm <- Y_1D[,-rm_loc]
#'   x_1Dnew <- as.matrix(seq(-5, 5, length = 100))
#'   cv_1D <- spatpca(x = x_1Drm, Y = Y_1Drm, tau2 = 1:100, x_new = x_1Dnew)
#'   plot(x_1Dnew, cv_1D$eigenfn, type = "l", main = "eigenfunction")
#'   plot(cv_1D$Yhat[, 50], xlab = "n", ylab = "Yhat", type = "l", 
#'        main = paste("prediction at x = ", x_1Dnew[50]))
#'   ###### 2D: Daily 8-hour ozone averages for sites in the Midwest (USA)
#'   library(fields)
#'   library(pracma)
#'   data(ozone2)
#'   x <- ozone2$lon.lat
#'   Y <- ozone2$y
#'   date <- as.Date(ozone2$date, format = "%y%m%d")
#'   rmna <- !colSums(is.na(Y))
#'   YY <- matrix(Y[, rmna], nrow = nrow(Y))
#'   YY <- detrend(YY, "linear")
#'   xx <- x[rmna, ]
#'   cv <- spatpca(x = xx, Y = YY)
#'   quilt.plot(xx, cv$eigenfn[,1])
#'   map("state", xlim = range(xx[, 1]), ylim = range(xx[, 2]), add = T)
#'   map.text("state", xlim = range(xx[, 1]), ylim = range(xx[, 2]), cex = 2, add = T)
#'   plot(date, YY %*% cv$eigenfn[,1], type = "l", ylab = "1st Principal Component")
#'   ### new loactions
#'   new_p = 200
#'   x_lon <- seq(min(xx[, 1]), max(xx[, 1]), length = new_p)
#'   x_lat <- seq(min(xx[, 2]), max(xx[, 2]), length = new_p)
#'   xx_new <- as.matrix(expand.grid(x = x_lon, y = x_lat))
#'   eof <- spatpca(x = xx, Y = YY, K = cv$Khat, tau1 = cv$stau1, tau2 = cv$stau2, x_new = xx_new)
#'   quilt.plot(xx_new, eof$eigenfn[,1], nx = new_p, ny = new_p, xlab = "lon.", ylab = "lat.")
#'   map("state", xlim = range(x_lon), ylim = range(x_lat), add = T)
#'   map.text("state", xlim = range(x_lon), ylim = range(x_lat), cex = 2, add = T)
#'   ###### 3D: regular locations
#'   x <- y <- z <- as.matrix(seq(-5, 5, length = 10))
#'   d <- expand.grid(x, y, z)
#'   Phi_3D <- exp(-d[, 1]^2 - d[, 2]^2 - d[, 3]^2) / norm(exp(-d[, 1]^2 - d[, 2]^2 - d[, 3]^2), "F")
#'   Y_3D <- rnorm(n = 1000, sd = 3) %*% t(Phi) + matrix(rnorm(n = 100 * 10^3), 100, 10^3)
#'   cv_3D <- spatpca(x = d, Y = Y_3D, tau2 = seq(0, 1000, length = 10))
#'   library(plot3D)
#'   library(RColorBrewer)
#'   cols <- colorRampPalette(brewer.pal(9, "Blues"))(10)
#'   isosurf3D(x, y, z, colvar = array(cv_3D$eigenfn[, 1], c(p, p, p)),
#'             level= seq(min(cv_3D$eigenfn[, 1]), max(cv_3D$eigenfn[, 1]), length = 10),
#'             ticktype = "detailed",
#'             colkey = list(side = 1),
#'             col = cols)
#' }
spatmca <- function(x1,
                    x2,
                    Y1,
                    Y2,
                    M = 5,
                    K = NULL,
                    doesSelectRank = ifelse(is.null(K), TRUE, FALSE),
                    tau1u = NULL,
                    tau2u = NULL,
                    tau1v = NULL,
                    tau2v = NULL,
                    x1New = NULL,
                    x2New = NULL,
                    center = TRUE,
                    plot.cv = FALSE,
                    maxit = 100,
                    thr = 1e-04,
                    doesSelectAllTuningParameters = FALSE,
                    numCores = NULL) {
  set_cores(numCores)
  
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  
  if (nrow(x1) != ncol(Y1))
    stop("The number of rows of x1 should be equal to the number of columns of Y1.")
  if (nrow(x1) < 3 || nrow(x2) < 3)
    stop("Number of locations must be larger than 2.")
  if (ncol(x1) > 3 || ncol(x2) > 3)
    stop("Dimension of locations must be less 4.")
  if (nrow(Y1) != nrow(Y2))
    stop("The numbers of sample sizes of both data should be equal.")
  if (M >= max(nrow(Y1)))
    stop("Number of folds must be less than sample size.")
  
  if (center == TRUE) {
    Y1 = Y1 - apply(Y1 , 2, "mean")
    Y2 = Y2 - apply(Y2 , 2, "mean")
  }
  n = nrow(Y1)
  stra <- sample(rep(1:M, length.out = nrow(Y1)))
  
  tempegvl1 <- svd(Y1 / n)
  tempegvl2 <- svd(Y2 / n)
  dd <- t(Y1) %*% Y2 / n
  tempegvl3 <- svd(dd)
  egvl1 <- tempegvl1$d[1]
  egvl2 <- tempegvl2$d[1]
  egvl3 <- tempegvl3$d[1]
  
  if (is.null(tau2u) && is.null(tau2v)) {
    ntau2u <- ntau2v <- 11
    
    indexu <-
      sort(abs(tempegvl3$u[, 1]),
           decreasing = T,
           index.return = T)$ix
    nu1u <- indexu[2]
    nu2u <- indexu[ncol(Y1)]
    max.tau2u <- 2 * abs(dd[nu1u, ] %*% tempegvl3$v[, 1])[1]
    min.tau2u <- abs(dd[nu2u, ] %*% tempegvl3$v[, 1])[1]
    
    tau2u <-
      c(0, exp(seq(
        log(min.tau2u), log(max.tau2u), length = (ntau2u - 1)
      )))
    
    indexv <-
      sort(abs(tempegvl3$v[, 1]),
           decreasing = T,
           index.return = T)$ix
    nu1v <- indexv[2]
    nu2v <- indexv[ncol(Y2)]
    max.tau2v <- 2 * abs(t(dd)[nu1v, ] %*% tempegvl3$u[, 1])[1]
    min.tau2v <- abs(t(dd)[nu2v, ] %*% tempegvl3$u[, 1])[1]
    tau2v <-
      c(0, exp(seq(
        log(min.tau2v), log(max.tau2v), length = (ntau2v - 1)
      )))
    
  } else if (is.null(tau2u)) {
    ntau2u <- 11
    indexu <-
      sort(abs(tempegvl3$u[, 1]),
           decreasing = T,
           index.return = T)$ix
    nu1u <- indexu[2]
    nu2u <- indexu[ncol(Y1)]
    max.tau2u <- 2 * abs(dd[nu1u, ] %*% tempegvl3$v[, 1])[1]
    min.tau2u <- abs(dd[nu2u, ] %*% tempegvl3$v[, 1])[1]
    tau2u <-
      c(0, exp(seq(
        log(min.tau2u), log(max.tau2u), length = (ntau2u - 1)
      )))
    
    ntau2v <- length(tau2v)
  } else if (is.null(tau2v)) {
    ntau2v <- 11
    indexv <-
      sort(abs(tempegvl3$v[, 1]),
           decreasing = T,
           index.return = T)$ix
    nu1v <- indexv[2]
    nu2v <- indexv[ncol(Y2)]
    max.tau2v <-
      egvl3 * abs(t(dd)[nu1v, ] %*% tempegvl3$u[, 1])[1]
    min.tau2v <-
      egvl3 * abs(t(dd)[nu2v, ] %*% tempegvl3$u[, 1])[1]
    tau2v <-
      c(0, exp(seq(
        log(min.tau2v), log(max.tau2v), length = (ntau2v - 1)
      )))
    
    ntau2u <- length(tau2u)
  } else{
    ntau2u <- length(tau2u)
    ntau2v <- length(tau2v)
  }
  
  
  if (is.null(tau1u) && is.null(tau1v)) {
    ntau1u <- 11
    ntau1v <- 11
    max.tau1u <- egvl3 / egvl1 * sqrt(ncol(Y1) / nrow(Y1))[1]
    max.tau1v <- egvl3 / egvl2 * sqrt(ncol(Y2) / nrow(Y2))[1]
    tau1u <-
      c(0, exp(seq(
        log(max.tau1u / 1e3), log(max.tau1u), length = (ntau1u - 1)
      )))
    tau1v <-
      c(0, exp(seq(
        log(max.tau1v / 1e3), log(max.tau1v), length = (ntau1v - 1)
      )))
  } else if (is.null(tau1u)) {
    ntau1u <- 11
    max.tau1u <- egvl3 / egvl1 * sqrt(ncol(Y1) / nrow(Y1))[1]
    ntau1v <- length(tau1v)
  } else if (is.null(tau1v)) {
    ntau1v <- 11
    max.tau1v <- egvl3 / egvl2 * sqrt(ncol(Y2) / nrow(Y2))[1]
    ntau1u <- length(tau1u)
  } else{
    ntau1u <- length(tau1u)
    ntau1v <- length(tau1v)
  }
  
  if (M < 2 && (max(ntau1u, ntau2u, ntau1v, ntau2v) > 1)) {
    ntau1u <- 1
    ntau2u <- 1
    ntau1v <- 1
    ntau2v <- 1
    warning("Only produce the result based on the largest tau1 and largest tau2.")
  }
  
  if (ntau2u == 1 && tau2u > 0) {
    if (tau2u != 0)
      l2u <-
        c(0, exp(seq(log(tau2u / 1e3), log(tau2u), length = 10)))
    else
      l2u <- tau2u
  } else{
    l2u <- 1
  }
  if (ntau2v == 1 && tau2v > 0) {
    if (tau2v != 0)
      l2v <-
        c(0, exp(seq(log(tau2v / 1e3), log(tau2v), length = 10)))
    else
      l2v <- tau2u
  }
  else{
    l2v <- 1
  }
  if (doesSelectRank == TRUE) {
    if (doesSelectAllTuningParameters == FALSE)
      cvtempold <- spatmcacv_rcpp(x1,
                                  x2,
                                  Y1,
                                  Y2,
                                  M,
                                  1,
                                  tau1u,
                                  tau2u,
                                  tau1v,
                                  tau2v,
                                  stra,
                                  maxit,
                                  thr,
                                  l2u,
                                  l2v)
    else{
      warning("Computing time may be quite long")
      cvtempold <- spatmcacvall_rcpp(x1,
                                     x2,
                                     Y1,
                                     Y2,
                                     M,
                                     1,
                                     tau1u,
                                     tau2u,
                                     tau1v,
                                     tau2v,
                                     stra,
                                     maxit,
                                     thr,
                                     l2u,
                                     l2v)
    }
    for (k in 2:min(dim(Y1), dim(Y2))) {
      if (doesSelectAllTuningParameters == FALSE)
        cvtemp <- spatmcacv_rcpp(x1,
                                 x2,
                                 Y1,
                                 Y2,
                                 M,
                                 k,
                                 tau1u,
                                 tau2u,
                                 tau1v,
                                 tau2v,
                                 stra,
                                 maxit,
                                 thr,
                                 l2u,
                                 l2v)
      else{
        warning("Computing time may be quite long")
        cvtemp <- spatmcacvall_rcpp(x1,
                                    x2,
                                    Y1,
                                    Y2,
                                    M,
                                    k,
                                    tau1u,
                                    tau2u,
                                    tau1v,
                                    tau2v,
                                    stra,
                                    maxit,
                                    thr,
                                    l2u,
                                    l2v)
      }
      
      if (min(cvtempold$cv2) <= min(cvtemp$cv2) ||
          abs(min(cvtempold$cv2) - min(cvtemp$cv2)) <= 1e-8)
        break
      cvtempold <- cvtemp
    }
    Khat <- k - 1
  }
  else{
    if (doesSelectAllTuningParameters == FALSE)
      cvtempold <- spatmcacv_rcpp(x1,
                                  x2,
                                  Y1,
                                  Y2,
                                  M,
                                  K,
                                  tau1u,
                                  tau2u,
                                  tau1v,
                                  tau2v,
                                  stra,
                                  maxit,
                                  thr,
                                  l2u,
                                  l2v)
    else{
      warning("Computing time may be quite long")
      cvtempold <- spatmcacvall_rcpp(x1,
                                     x2,
                                     Y1,
                                     Y2,
                                     M,
                                     K,
                                     tau1u,
                                     tau2u,
                                     tau1v,
                                     tau2v,
                                     stra,
                                     maxit,
                                     thr,
                                     l2u,
                                     l2v)
    }
    Khat <- K
  }
  
  cvtau1u <- cvtempold$cvtau1u
  cvtau2u <- cvtempold$cvtau2u
  cvtau1v <- cvtempold$cvtau1v
  cvtau2v <- cvtempold$cvtau2v
  cv1 <- cvtempold$cv1
  cv2 <- cvtempold$cv2
  cvall <- cvtempold$cvall
  Uest <- cvtempold$Uest
  Vest <- cvtempold$Vest
  if (is.null(x1New)) {
    x1New = x1
    Uestfn <- Uest
  }
  else{
    x1New = as.matrix(x1New)
    Uestfn <- tpm2(x1New, x1, Uest)
  }
  if (is.null(x2New)) {
    x2New = x2
    Vestfn <- Vest
  }
  else{
    x2New = as.matrix(x2New)
    Vestfn <- tpm2(x2New, x2, Vest)
  }
  if (plot.cv == TRUE && !is.null(cv1)) {
    par(mfrow = c(2, 1))
    image.plot(tau1u, tau1v, cv1, main = "for tau1u and tau1v selection given tau2u and tau2v")
    image.plot(tau2u, tau2v, cv2, main = "for tau2u and tau2v selection given selected tau1u and tau2v")
  }
  Dest <- as.vector(cvtempold$Dest)
  crosscovfn <-
    Uestfn %*% diag(Dest, nrow = Khat, ncol = Khat) %*% t(Vestfn)
  obj.cv <-
    list(
      Uestfn = Uestfn,
      Vestfn = Vestfn,
      crosscov = crosscovfn,
      Dest = Dest,
      cv1 = cv1,
      cv2 = cv2,
      cvall = cvall,
      Khat = Khat,
      stau1u = cvtau1u,
      stau2u = cvtau2u,
      stau1v = cvtau1v,
      stau2v = cvtau2v,
      tau1u = tau1u,
      tau2u = tau2u,
      tau1v = tau1v,
      tau2v = tau2v
    )
  return(obj.cv)
}
