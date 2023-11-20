#' @title  Regularized spatial MCA
#'
#' @description Produce spatial coupled patterns at the designated locations according to the specified tuning parameters or the tuning parameters selected by M-fold cross-validation.
#'
#' @param x1 Location matrix (\eqn{p \times d}) corresponding to Y1. Each row is a location. \eqn{d=1,2} is the dimension of locations.
#' @param x2 Location matrix (\eqn{q \times d}) corresponding to Y2. Each row is a location.
#' @param Y1 Data matrix (\eqn{n \times p}) of the first variable stores the values at \eqn{p} locations with sample size \eqn{n}.
#' @param Y2 Data matrix (\eqn{n \times q}) of the second variable stores the values at \eqn{q} locations with sample size \eqn{n}.
#' @param M Optional number of folds; default is 5.
#' @param K Optional user-supplied number of coupled patterns; default is NULL. If K is NULL or is_K_selected is TRUE, K is selected automatically.
#' @param is_K_selected If TRUE, K is selected automatically; otherwise, is_K_selected is set to be user-supplied K. Default depends on user-supplied K.
#' @param tau1u Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence corresponding to Y1. If NULL, 10 tau1u values in a range are used.
#' @param tau2u Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence corresponding to Y1. If NULL, 10 tau2u values in a range are used.
#' @param tau1v Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence corresponding to Y2. If NULL, 10 tau1v values in a range are used.
#' @param tau2v Optional user-supplied numeric vector of a nonnegative smoothness parameter sequence corresponding to Y2. If NULL, 10 tau2v values in a range are used.
#' @param x1New New location matrix corresponding to Y1. If NULL, it is x1.
#' @param x2New New location matrix corresponding to Y2. If NULL, it is x2.
#' @param center If TRUE, center the columns of Y. Default is FALSE.
#' @param maxit Maximum number of iterations. Default value is 100.
#' @param thr Threshold for convergence. Default value is \eqn{10^{-4}}.
#' @param are_all_tuning_parameters_selected If TRUE, the K-fold CV performs to select 4 tuning parameters simultaneously. Default value is FALSE.
#' @param num_cores Number of cores used to parallel computing. Default value is NULL (See `RcppParallel::defaultNumThreads()`)
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
#' \item{cv1}{cv scores for tau1u and tau1v when are_all_tuning_parameters_selected is FALSE.}
#' \item{cv2}{cv scores for tau2u and tau2v when are_all_tuning_parameters_selected is FALSE.}
#' \item{cvall}{cv scores for tau1u, tau2u, tau1v and tau2v when are_all_tuning_parameters_selected is TRUE.}
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
#' originalPar <- par(no.readonly = TRUE)
#' # The following examples only use two threads for parallel computing.
#' ## 1D: regular locations
#' p <- q <- 10
#' n <- 100
#' x1 <- matrix(seq(-7, 7, length = p), nrow = p, ncol = 1)
#' x2 <- matrix(seq(-7, 7, length = q), nrow = q, ncol = 1)
#' u <- exp(-x1^2) / norm(exp(-x1^2), "F")
#' v <- exp(-(x2 - 2)^2) / norm(exp(-(x2 - 2)^2), "F")
#' Sigma <- array(0, c(p + q, p + q))
#' Sigma[1:p, 1:p] <- diag(p)
#' Sigma[(p + 1):(p + q), (p + 1):(p + q)] <- diag(p)
#' Sigma[1:p, (p + 1):(p + q)] <- u %*% t(v)
#' Sigma[(p + 1):(p + q), 1:p] <- t(Sigma[1:p, (p + 1):(p + q)])
#' noise <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = 0.001 * diag(p + q))
#' Y <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = Sigma) + noise
#' Y1 <- Y[, 1:p]
#' Y2 <- Y[, -(1:p)]
#' cv1 <- spatmca(x1, x2, Y1, Y2, num_cores = 2)
#'
#' par(mfrow = c(2, 1))
#' plot(x1, cv1$Uestfn[, 1], type='l', main = "1st pattern for Y1")
#' plot(x1, cv1$Vestfn[, 1], type='l', main = "1st pattern for Y2")
#' ## Avoid changing the global enviroment
#' par(originalPar)
#'
#' \donttest{
#' # The following examples will be executed more than 5 secs or including other libraries.
#' ## 1D: artificial irregular locations
#' rmLoc1 <- sample(1:p, 3)
#' rmLoc2 <- sample(1:q, 4)
#' x1Rm <- x1[-rmLoc1]
#' x2Rm <- x2[-rmLoc2]
#' Y1Rm <- Y1[, -rmLoc1]
#' Y2Rm <- Y2[, -rmLoc2]
#' x1New <- as.matrix(seq(-7, 7, length = 100))
#' x2New <- as.matrix(seq(-7, 7, length = 50))
#' cv2 <- spatmca(x1 = x1Rm,
#'                x2 = x2Rm,
#'                Y1 = Y1Rm,
#'                Y2 = Y2Rm,
#'                x1New = x1New,
#'                x2New = x2New)
#' par(mfrow = c(2, 1))
#' plot(x1New, cv2$Uestfn[,1], type='l', main = "1st pattern for Y1")
#' plot(x2New, cv2$Vestfn[,1], type='l', main = "1st pattern for Y2")
#' par(originalPar)
#'
#' ## 2D real data
#' ##  Daily 8-hour ozone averages and maximum temperature obtained from 28 monitoring
#' ##  sites of NewYork, USA. It is of interest to see the relationship between the ozone
#' ##  and the temperature through the coupled patterns.
#'
#' library(spTimer)
#' library(pracma)
#' library(fields)
#' library(maps)
#' data(NYdata)
#' NYsite <- unique(cbind(NYdata[, 1:3]))
#' date <- as.POSIXct(seq(as.Date("2006-07-01"), as.Date("2006-08-31"), by = 1))
#' cMAXTMP<- matrix(NYdata[,8], 62, 28)
#' oz <- matrix(NYdata[,7], 62, 28)
#' rmNa <- !colSums(is.na(oz))
#' temp <- detrend(matrix(cMAXTMP[, rmNa], nrow = nrow(cMAXTMP)), "linear")
#' ozone <- detrend(matrix(oz[, rmNa], nrow = nrow(oz)), "linear")
#' x1 <- NYsite[rmNa, 2:3]
#' cv <- spatmca(x1, x1, temp, ozone)
#' par(mfrow = c(2, 1))
#' quilt.plot(x1, cv$Uestfn[, 1],
#'            xlab = "longitude",
#'            ylab = "latitude",
#'            main = "1st spatial pattern for temperature")
#' map(database = "state", regions = "new york", add = TRUE)
#' quilt.plot(x1, cv$Vestfn[, 1],
#'            xlab = "longitude",
#'            ylab = "latitude",
#'            main = "1st spatial pattern for ozone")
#' map(database = "state", regions = "new york", add = TRUE)
#' par(originalPar)
#'
#' ### Time series for the coupled patterns
#' tstemp <- temp %*% cv$Uestfn[,1]
#' tsozone <- ozone %*% cv$Vestfn[,1]
#' corr <- cor(tstemp, tsozone)
#' plot(date, tstemp / sd(tstemp), type='l', main = "Time series", ylab = "", xlab = "month")
#' lines(date, tsozone/sd(tsozone),col=2)
#' legend("bottomleft", c("Temperature (standardized)", "Ozone (standardized)"), col = 1:2, lty = 1:1)
#' mtext(paste("Pearson's correlation = ", round(corr, 3)), 3)
#'
#  ### New locations
#' newP <- 50
#' xLon <- seq(-80, -72, length = newP)
#' xLat <- seq(41, 45, length = newP)
#' xxNew <- as.matrix(expand.grid(x = xLon, y = xLat))
#' cvNew <- spatmca(x1 = x1,
#'                  x2 = x1,
#'                  Y1 = temp,
#'                  Y2 = ozone,
#'                  K = cv$Khat,
#'                  tau1u = cv$stau1u,
#'                  tau1v = cv$stau1v,
#'                  tau2u = cv$stau2u,
#'                  tau2v = cv$stau2v,
#'                  x1New = xxNew,
#'                  x2New = xxNew)
#' par(mfrow = c(2, 1))
#' quilt.plot(xxNew, cvNew$Uestfn[, 1],
#'            nx = newP,
#'            ny = newP,
#'            xlab = "longitude",
#'            ylab = "latitude",
#'            main = "1st spatial pattern for temperature")
#' map(database = "county", regions = "new york", add = TRUE)
#' map.text("state", regions = "new york", cex = 2, add = TRUE)
#' quilt.plot(xxNew, cvNew$Vestfn[, 1],
#'            nx = newP,
#'            ny = newP,
#'            xlab = "longitude",
#'            ylab = "latitude",
#'            main = "2nd spatial pattern for ozone")
#' map(database = "county", regions = "new york", add = TRUE)
#' map.text("state", regions = "new york", cex = 2, add = TRUE)
#' par(originalPar)
#'
#' ## 3D: regular locations
#' n <- 200
#' x <- y <- z <- as.matrix(seq(-7, 7, length = 8))
#' d <- expand.grid(x, y, z)
#' u3D <- v3D <- exp(-d[, 1]^2 - d[, 2]^2 -d[, 3]^2)
#' p <- q <- 8^3
#' Sigma3D <- array(0, c(p + q, p + q))
#' Sigma3D[1:p, 1:p] <- diag(p)
#' Sigma3D[(p + 1):(p + q), (p + 1):(p + q)] <- diag(p)
#' Sigma3D[1:p, (p + 1):(p + q)] <- u3D %*% t(v3D)
#' Sigma3D[(p + 1):(p + q), 1:p] <- t(Sigma3D[1:p, (p + 1):(p + q)])
#'
#' noise3D <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = 0.001 * diag(p + q))
#' Y3D <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = Sigma3D) + noise3D
#' Y13D <- Y3D[, 1:p]
#' Y23D <- Y3D[, -(1:p)]
#' cv3D <- spatmca(d, d, Y13D, Y23D)
#'
#' library(plot3D)
#' library(RColorBrewer)
#' cols <- colorRampPalette(brewer.pal(9, 'Blues'))(10)
#' isosurf3D(x, y, z,
#'           colvar = array(cv3D$Uestfn[, 1], c(8, 8, 8)),
#'           level = seq(min(cv3D$Uestfn[, 1]), max(cv3D$Uestfn[, 1]), length = 10),
#'           ticktype = "detailed",
#'           colkey = list(side = 1),
#'           col = cols,
#'           main = "1st estimated pattern for Y1")
#'
#' isosurf3D(x, y, z,
#'           colvar = array(cv3D$Vestfn[, 1], c(8, 8, 8)),
#'           level = seq(min(cv3D$Vestfn[, 1]), max(cv3D$Vestfn[,1]), length = 10),
#'           ticktype = "detailed",
#'           colkey = list(side = 1),
#'           col = cols,
#'           main = "1st estimated pattern for Y2")
#' }
spatmca <- function(x1,
                    x2,
                    Y1,
                    Y2,
                    M = 5,
                    K = NULL,
                    is_K_selected = ifelse(is.null(K), TRUE, FALSE),
                    tau1u = NULL,
                    tau2u = NULL,
                    tau1v = NULL,
                    tau2v = NULL,
                    x1New = NULL,
                    x2New = NULL,
                    center = TRUE,
                    maxit = 100,
                    thr = 1e-04,
                    are_all_tuning_parameters_selected = FALSE,
                    num_cores = NULL) {
  call <- match.call()
  setCores(num_cores)
  
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  checkInputData(x1, x2, Y1, Y2, M)
  
  Y1 <- detrend(Y1, is_K_selected)
  Y2 <- detrend(Y2, is_K_selected)
  
  n <- nrow(Y1)
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
           decreasing = TRUE,
           index.return = TRUE)$ix
    nu1u <- indexu[2]
    nu2u <- indexu[ncol(Y1)]
    max.tau2u <- 2 * abs(dd[nu1u, ] %*% tempegvl3$v[, 1])[1]
    min.tau2u <- abs(dd[nu2u, ] %*% tempegvl3$v[, 1])[1]
    
    tau2u <-
      c(0, exp(seq(
        log(min.tau2u), log(max.tau2u),
        length = (ntau2u - 1)
      )))
    
    indexv <-
      sort(abs(tempegvl3$v[, 1]),
           decreasing = TRUE,
           index.return = TRUE)$ix
    nu1v <- indexv[2]
    nu2v <- indexv[ncol(Y2)]
    max.tau2v <- 2 * abs(t(dd)[nu1v, ] %*% tempegvl3$u[, 1])[1]
    min.tau2v <- abs(t(dd)[nu2v, ] %*% tempegvl3$u[, 1])[1]
    tau2v <-
      c(0, exp(seq(
        log(min.tau2v), log(max.tau2v),
        length = (ntau2v - 1)
      )))
  } else if (is.null(tau2u)) {
    ntau2u <- 11
    indexu <-
      sort(abs(tempegvl3$u[, 1]),
           decreasing = TRUE,
           index.return = TRUE)$ix
    nu1u <- indexu[2]
    nu2u <- indexu[ncol(Y1)]
    max.tau2u <- 2 * abs(dd[nu1u, ] %*% tempegvl3$v[, 1])[1]
    min.tau2u <- abs(dd[nu2u, ] %*% tempegvl3$v[, 1])[1]
    tau2u <-
      c(0, exp(seq(
        log(min.tau2u), log(max.tau2u),
        length = (ntau2u - 1)
      )))
    
    ntau2v <- length(tau2v)
  } else if (is.null(tau2v)) {
    ntau2v <- 11
    indexv <-
      sort(abs(tempegvl3$v[, 1]),
           decreasing = TRUE,
           index.return = TRUE)$ix
    nu1v <- indexv[2]
    nu2v <- indexv[ncol(Y2)]
    max.tau2v <-
      egvl3 * abs(t(dd)[nu1v, ] %*% tempegvl3$u[, 1])[1]
    min.tau2v <-
      egvl3 * abs(t(dd)[nu2v, ] %*% tempegvl3$u[, 1])[1]
    tau2v <-
      c(0, exp(seq(
        log(min.tau2v), log(max.tau2v),
        length = (ntau2v - 1)
      )))
    
    ntau2u <- length(tau2u)
  } else {
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
        log(max.tau1u / 1e3), log(max.tau1u),
        length = (ntau1u - 1)
      )))
    tau1v <-
      c(0, exp(seq(
        log(max.tau1v / 1e3), log(max.tau1v),
        length = (ntau1v - 1)
      )))
  } else if (is.null(tau1u)) {
    ntau1u <- 11
    max.tau1u <- egvl3 / egvl1 * sqrt(ncol(Y1) / nrow(Y1))[1]
    ntau1v <- length(tau1v)
    tau1u <-
      c(0, exp(seq(
        log(max.tau1u / 1e3), log(max.tau1u),
        length = (ntau1u - 1)
      )))
  } else if (is.null(tau1v)) {
    ntau1v <- 11
    max.tau1v <- egvl3 / egvl2 * sqrt(ncol(Y2) / nrow(Y2))[1]
    ntau1u <- length(tau1u)
    tau1v <-
      c(0, exp(seq(
        log(max.tau1v / 1e3), log(max.tau1v),
        length = (ntau1v - 1)
      )))
  } else {
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
  
  if (ntau2u == 1 && length(tau2u) == 1) {
    if (tau2u != 0) {
      l2u <-
        c(0, exp(seq(
          log(tau2u / 1e3), log(tau2u), length = 10
        )))
    } else {
      l2u <- tau2u
    }
  } else {
    l2u <- 1
  }
  if (ntau2v == 1 && length(tau2v) == 1) {
    if (tau2v != 0) {
      l2v <-
        c(0, exp(seq(
          log(tau2v / 1e3), log(tau2v), length = 10
        )))
    } else{
      l2v <- tau2u
    }
  } else {
    l2v <- 1
  }
  if (is_K_selected) {
    if (are_all_tuning_parameters_selected == FALSE) {
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
    } else {
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
      if (are_all_tuning_parameters_selected == FALSE) {
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
      } else {
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
          abs(min(cvtempold$cv2) - min(cvtemp$cv2)) <= 1e-8) {
        break
      }
      cvtempold <- cvtemp
    }
    Khat <- k - 1
  }
  else {
    if (are_all_tuning_parameters_selected == FALSE) {
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
    } else {
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
    x1New <- x1
    Uestfn <- Uest
  }
  else {
    x1New <- as.matrix(x1New)
    Uestfn <- tpm2(x1New, x1, Uest)
  }
  if (is.null(x2New)) {
    x2New <- x2
    Vestfn <- Vest
  }
  else {
    x2New <- as.matrix(x2New)
    Vestfn <- tpm2(x2New, x2, Vest)
  }
  Dest <- as.vector(cvtempold$Dest)
  crosscovfn <-
    Uestfn %*% diag(Dest, nrow = Khat, ncol = Khat) %*% t(Vestfn)
  obj.cv <-
    list(
      call = call,
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
  class(obj.cv) <- "spatmca"
  return(obj.cv)
}


#'
#' @title  Display the cross-validation results
#'
#' @description Display the M-fold cross-validation results
#'
#' @param x An spatmca class object for `plot` method
#' @param ... Not used directly
#' @return `NULL`
#' @seealso \link{spatmca}
#'
#' @export
#' @method plot spatmca
#' @examples
#' p <- q <- 5
#' n <- 50
#' x1 <- matrix(seq(-7, 7, length = p), nrow = p, ncol = 1)
#' x2 <- matrix(seq(-7, 7, length = q), nrow = q, ncol = 1)
#' u <- exp(-x1^2) / norm(exp(-x1^2), "F")
#' v <- exp(-(x2 - 2)^2) / norm(exp(-(x2 - 2)^2), "F")
#' Sigma <- array(0, c(p + q, p + q))
#' Sigma[1:p, 1:p] <- diag(p)
#' Sigma[(p + 1):(p + q), (p + 1):(p + q)] <- diag(p)
#' Sigma[1:p, (p + 1):(p + q)] <- u %*% t(v)
#' Sigma[(p + 1):(p + q), 1:p] <- t(Sigma[1:p, (p + 1):(p + q)])
#' noise <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = 0.001 * diag(p + q))
#' Y <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = Sigma) + noise
#' Y1 <- Y[, 1:p]
#' Y2 <- Y[, -(1:p)]
#' cv_1D <- spatmca(x1, x2, Y1, Y2, num_cores = 2)
#' plot(cv_1D)
#
plot.spatmca <- function(x, ...) {
  if (!inherits(x, "spatmca")) {
    stop("Invalid object! Please enter a `spatmca` object")
  }
  
  cv_data <- result <- list()
  variate_names <- c("First Variate", "Second Variate")
  
  cv_data[[variate_names[1]]] <-
    expand.grid(u = x$tau1u, v = x$tau1v)
  cv_data[[variate_names[1]]]$cv <- as.vector(x$cv1)
  cv_data[[variate_names[2]]] <-
    expand.grid(u = x$tau2u, v = x$tau2v)
  cv_data[[variate_names[2]]]$cv <- as.vector(x$cv2)
  
  for (variate in variate_names) {
    result[[variate]] <- plot_cv_field(cv_data[[variate]], variate)
  }
  plot_sequentially(result)
}