spatmca <- function(x1, x2, Y1, Y2, M = 5, K = NULL, K.select = ifelse(is.null(K),TRUE,FALSE), tau1u = NULL, tau2u = NULL,
                    tau1v = NULL, tau2v = NULL, x1_new = NULL, x2_new = NULL, center = TRUE, plot.cv = FALSE, maxit = 100, thr = 1e-04, fullselect = FALSE){
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  
  if(nrow(x1) != ncol(Y1))
    stop("The number of rows of x1 should be equal to the number of columns of Y1.")
  if (nrow(x1) < 3 ||nrow(x2) < 3)
    stop("Number of locations must be larger than 2.")
  if (ncol(x1) > 3 || ncol(x2) > 3)
    stop("Dimension of locations must be less 4.")
  if(nrow(Y1) != nrow(Y2))
    stop("The numbers of sample sizes of both data should be equal.")
  if(M >= max(nrow(Y1)))
    stop("Number of folds must be less than sample size.")
  
  if(center == TRUE){
    Y1 = Y1 - apply(Y1 , 2, "mean")
    Y2 = Y2 - apply(Y2 , 2, "mean")
  }
  n = nrow(Y1)
  stra <- sample(rep(1:M, length.out = nrow(Y1)))
  
  tempegvl1 <- svd(Y1/n)
  tempegvl2 <- svd(Y2/n)
  dd <- t(Y1)%*%Y2/n
  tempegvl3 <- svd(dd)
  egvl1 <- tempegvl1$d[1]
  egvl2 <- tempegvl2$d[1]
  egvl3 <- tempegvl3$d[1]
  
  if(is.null(tau2u)&&is.null(tau2v)){
    ntau2u <- ntau2v <- 11
    
    indexu <- sort(abs(tempegvl3$u[,1]), decreasing=T, index.return=T)$ix
    nu1u <- indexu[2]
    nu2u <- indexu[ncol(Y1)]
    max.tau2u <- 2*abs(dd[nu1u,]%*%tempegvl3$v[,1])[1]
    min.tau2u <- abs(dd[nu2u,]%*%tempegvl3$v[,1])[1]
  
    tau2u <- c(0,exp(seq(log(min.tau2u), log(max.tau2u), length = (ntau2u-1)))) 
    
    indexv <- sort(abs(tempegvl3$v[,1]),decreasing=T,index.return=T)$ix
    nu1v <- indexv[2]
    nu2v <- indexv[ncol(Y2)]
    max.tau2v <- 2*abs(t(dd)[nu1v,]%*%tempegvl3$u[,1])[1]
    min.tau2v <- abs(t(dd)[nu2v,]%*%tempegvl3$u[,1])[1]
    tau2v <- c(0,exp(seq(log(min.tau2v), log(max.tau2v), length = (ntau2v-1)))) 
    
  }else if (is.null(tau2u)){ 
    ntau2u <- 11
    indexu <- sort(abs(tempegvl3$u[,1]),decreasing=T,index.return=T)$ix
    nu1u <- indexu[2]
    nu2u <- indexu[ncol(Y1)]
    max.tau2u <- 2*abs(dd[nu1u,]%*%tempegvl3$v[,1])
    min.tau2u <- abs(dd[nu2u,]%*%tempegvl3$v[,1])
    tau2u <- c(0,exp(seq(log(min.tau2u), log(max.tau2u), length = (ntau2u-1)))) 
    
    ntau2v <- length(tau2v)
  }else if (is.null(tau2v)){ 
    ntau2v <- 11
    indexv <- sort(abs(tempegvl3$v[,1]),decreasing=T,index.return=T)$ix
    nu1v <- indexv[2]
    nu2v <- indexv[ncol(Y2)]
    max.tau2v <- egvl3*abs(t(dd)[nu1v,]%*%tempegvl3$u[,1])
    min.tau2v <- egvl3*abs(t(dd)[nu2v,]%*%tempegvl3$u[,1])
    tau2v <- c(0,exp(seq(log(min.tau2v), log(max.tau2v), length = (ntau2v-1)))) 
    
    ntau2u <- length(tau2u)
  }else{
    ntau2u <- length(tau2u)
    ntau2v <- length(tau2v)
  }
  
  
  if(is.null(tau1u) && is.null(tau1v)) {
    ntau1u <- 11
    ntau1v <- 11
    max.tau1u <- egvl3/egvl1*sqrt(ncol(Y1)/nrow(Y1))
    max.tau1v <- egvl3/egvl2*sqrt(ncol(Y2)/nrow(Y2))
    tau1u <- c(0,exp(seq(log(max.tau1u/1e3), log(max.tau1u), length = (ntau1u-1)))) 
    tau1v <- c(0,exp(seq(log(max.tau1v/1e3), log(max.tau1v), length = (ntau1v-1)))) 
  }else if (is.null(tau1u)){
    ntau1u <- 11
    max.tau1u <- egvl3/egvl1*sqrt(ncol(Y1)/nrow(Y1))
    ntau1v <- length(tau1v)
  }else if (is.null(tau1v)){
    ntau1v <- 11
    max.tau1v <- egvl3/egvl2*sqrt(ncol(Y2)/nrow(Y2))
    ntau1u <- length(tau1u)
  }else{
    ntau1u <- length(tau1u)
    ntau1v <- length(tau1v)
  }
  
  if(M < 2 && (max(ntau1u, ntau2u, ntau1v,ntau2v) > 1)) {
    ntau1u <- 1
    ntau2u <- 1
    ntau1v <- 1
    ntau2v <- 1
    warning("Only produce the result based on the largest tau1 and largest tau2.")
  }  
  
  if(ntau2u == 1 && tau2u > 0){
    if(tau2u !=0)
      l2u <- c(0,exp(seq(log(tau2u/1e3), log(tau2u), length = 10)))
    else
      l2u <- tau2u
  }else{
    l2u <- 1
  }
  if(ntau2v == 1 && tau2v > 0){
    if(tau2v !=0)
      l2v <- c(0,exp(seq(log(tau2v/1e3), log(tau2v), length = 10)))
    else
      l2v <- tau2u
  }
  else{
    l2v <- 1
  }
  if(K.select == TRUE){
    if(fullselect == FALSE)
      cvtempold <- spatmcacv_rcpp(x1, x2, Y1, Y2, M, 1, tau1u, tau2u, tau1v, tau2v,  stra, maxit, thr, l2u, l2v)
    else{
      warning("Computing time may be quite long")
      cvtempold <- spatmcacvall_rcpp(x1, x2, Y1, Y2, M, 1, tau1u, tau2u, tau1v, tau2v,  stra, maxit, thr, l2u, l2v)
    }
    for(k in 2:min(dim(Y1),dim(Y2))){
      if(fullselect == FALSE)
        cvtemp <- spatmcacv_rcpp(x1, x2, Y1, Y2, M, k, tau1u, tau2u, tau1v, tau2v,  stra, maxit, thr, l2u, l2v)
      else{
        warning("Computing time may be quite long")
        cvtemp <- spatmcacvall_rcpp(x1, x2, Y1, Y2, M, k, tau1u, tau2u, tau1v, tau2v,  stra, maxit, thr, l2u, l2v)
      }
        
      if(min(cvtempold$cv2)<= min(cvtemp$cv2)||abs(min(cvtempold$cv2) - min(cvtemp$cv2))<=1e-8)
        break
      cvtempold <- cvtemp
    }
    Khat <- k-1
  }
  else{
    if(fullselect == FALSE)
      cvtempold <- spatmcacv_rcpp(x1, x2, Y1, Y2, M, K, tau1u, tau2u, tau1v, tau2v,  stra, maxit, thr, l2u, l2v)
    else{
      warning("Computing time may be quite long")
      cvtempold <- spatmcacvall_rcpp(x1, x2, Y1, Y2, M, K, tau1u, tau2u, tau1v, tau2v,  stra, maxit, thr, l2u, l2v)
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
  if(is.null(x1_new)){
    x1_new = x1
    Uestfn <- Uest
  }
  else{
    x1_new = as.matrix(x1_new)
    Uestfn <- tpm2(x1_new, x1, Uest)
  }
  if(is.null(x2_new)){
    x2_new = x2
    Vestfn <- Vest
  }
  else{
    x2_new = as.matrix(x2_new)
    Vestfn <- tpm2(x2_new, x2, Vest)
  }
  if(plot.cv == TRUE && !is.null(cv1)){
    par(mfrow=c(2,1))
    image.plot(tau1u,tau1v,cv1, main="for tau1u and tau1v selection given tau2u and tau2v")
    image.plot(tau2u, tau2v, cv2, main="for tau2u and tau2v selection given selected tau1u and tau2v")
  }
  Dest <- as.vector(cvtempold$Dest)
  crosscovfn <- Uestfn%*%diag(Dest, nrow = Khat, ncol = Khat)%*%t(Vestfn)
  obj.cv <- list(Uestfn = Uestfn, Vestfn = Vestfn, crosscov = crosscovfn, Dest = Dest, cv1 = cv1, cv2 = cv2, cvall = cvall, Khat = Khat,
                 stau1u = cvtau1u, stau2u = cvtau2u,stau1v = cvtau1v, stau2v = cvtau2v,
                 tau1u = tau1u, tau2u = tau2u, tau1v = tau1v, tau2v = tau2v)
  return(obj.cv)
}
