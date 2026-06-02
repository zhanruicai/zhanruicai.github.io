# --------------------------------------------
# Load required packages
# --------------------------------------------
library(geigen)
library(VGAM)
library(Matrix)
library(foreach)
library(parallel)
library(doParallel)
library(MASS)

# Load all functions and variables defined in functions.R into the current R session
source("functions.R")

# --------------------------------------------
# Function: dp_sparse_sir
# Purpose: Perform sparse differentially private SIR with two DP refinement methods
# Inputs:
#   X         - predictors
#   Y         - response
#   H         - number of slices
#   K         - number of dimensions to estimate
#   s         - sparsity level (number of active features)
#   Cn        - penalty parameter for model selection
#   TT1, TT2  - number of iterations for DP Rayleigh Flow and projected DP-SIR
#   stepp     - learning rate for projected DP-SIR
#   lambda    - regularization strength
# Outputs:
#   A vector of estimation errors from
# --------------------------------------------

dp_sparse_sir<-function(X, Y, H, K, s,Cn ,TT1,TT2,stepp,lambda){
  # Step 1: Slice the response using differentially private histogram
  n <- dim(X)[1]
  p <- dim(X)[2]
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  
  # Step 2: Estimate covariance matrix and SIR matrix
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  
  # Step 3: Oracle SIR (on top-s active features)
  sir.estimator <- geigen(A = Martix.M[1:s,1:s], B = Sigma.X[1:s,1:s], symmetric = T)
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.sir <- norm(proj_mat(beta1[1:s,])-proj_mat(sir.vector),type = "F")
  
  # Step 4: Initial estimate with noisy eigen-decomposition
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = s,epsilon = 3.5, delta = delta,k=K)
  sigma.e1 <- 2*s/2/n*sqrt(2*log(2.5/delta))
  # Add noise to matrices
  E1 <- rnorm(n = (s+1)*s/2, mean = 0, sd = sigma.e1)
  E1_noise <- matrix(data = 0, nrow = s, ncol = s)
  E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
  E1_noise <- E1_noise+t(E1_noise)
  diag(E1_noise)<- diag(E1_noise)/2
  sigma.e2 <- 2*s/2/n*sqrt(2*log(2.5/delta))
  E2 <- rnorm(n = (s+1)*s/2, mean = 0, sd = sigma.e2)
  E2_noise <- matrix(data = 0, nrow = s, ncol = s)
  E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E2
  E2_noise <- E2_noise+t(E2_noise)
  diag(E2_noise)<- diag(E2_noise)/2
  Sigma.X.noise.PD <- as.matrix(nearPD(Sigma.X[pindex,pindex] + E1_noise)$mat)
  ini.sir.estimator <- geigen(A = Martix.M[pindex,pindex] + E2_noise, B = Sigma.X.noise.PD, symmetric = T)
  
  order.lambda <- ini.sir.estimator$values[order(ini.sir.estimator$values,decreasing = T)]
  K = dp_bic(order.lambda = order.lambda, Cn = Cn, HH = length(unique(tildey)))
  
  ini.sir.number <- order(ini.sir.estimator$values,decreasing = T)[1:K]
  ini.sir.vector <- ini.sir.estimator$vectors[,ini.sir.number]
  ini.sir.vector.complete <-matrix(data = 0, nrow = p, ncol = K)
  ini.sir.vector.complete[pindex,] <- ini.sir.vector
  error.ini <- norm(proj_mat(beta1)-proj_mat((ini.sir.vector.complete)),type = "F")
  
  # Step 5: DP Rayleigh Flow refinement
  
  epsilon <- 1
  delta <- 1/n^{1.1}
  TT <- TT1
  if (K==1) {dp.rf <- t(ini.sir.vector.complete)} else {dp.rf <- (ini.sir.vector.complete[,1]);dim(dp.rf)<-c(1,p);}
  dim(dp.rf) <- c(1,p)
  for (T.iter in 1:TT) {
    sigma.e1 <- 2*s/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
    sigma.e2 <- 2*6.01*s/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
    E1 <- rnorm(n = (p+1)*p/2, mean = 0, sd = sigma.e1)
    E1_noise <- matrix(data = 0, nrow = p, ncol = p)
    E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
    E1_noise <- E1_noise+t(E1_noise)
    diag(E1_noise)<- diag(E1_noise)/2
    E2 <- rnorm(n = (p+1)*p/2, mean = 0, sd = sigma.e2)
    E2_noise <- matrix(data = 0, nrow = p, ncol = p)
    E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E2
    E2_noise <- E2_noise+t(E2_noise)
    diag(E2_noise)<- diag(E2_noise)/2
    Sigma.X.noise.PD <- Sigma.X + E1_noise
    rhot <- as.numeric((dp.rf)%*%(Martix.M + E2_noise)%*%t(dp.rf)/((dp.rf)%*%as.matrix( Sigma.X.noise.PD)%*%t(dp.rf)))
    Ct <- diag(p)+((Martix.M + E2_noise)-rhot* Sigma.X.noise.PD)/rhot/20 ##step size:1/2
    dp.rf <- dp.rf%*%Ct/norm(dp.rf%*%Ct,type = "F")
    dp.rf[-order(abs(dp.rf),decreasing = T)[1:s]]<-0
  }
  error.ray <- norm(proj_mat(beta1)-proj_mat(t(dp.rf)),type = "F")
  
  # Step 6: Projected DP-SIR with sparsity constraint
  if (K==1) {dp.sir <- ini.sir.vector.complete*sqrt(max(ini.sir.estimator$values)/lambda+1)}
  if (K>1) {dp.sir <- ini.sir.vector.complete%*%diag(sqrt(ini.sir.estimator$values[order(ini.sir.estimator$values,decreasing = T)[1:K]]/lambda+1))}
  dim(dp.sir)<-c(p,K)
  TT <- TT2
  for (T.iter in 1:TT) {
    sigma.e1 <- (6.01*4+lambda*2*4)/n*sqrt(K)/stepp
    gd <- -Martix.M%*%dp.sir+lambda*Sigma.X%*%dp.sir%*%(t(dp.sir)%*%Sigma.X%*%dp.sir-diag(K))
    gd <- gd
    E1 <- rnorm(n = s*K, mean = 0, sd = 6/5*sigma.e1*sqrt(s)*sqrt(2*log(1.25*TT*6/5/delta))/epsilon*TT)
    dp.sir <- dp.sir - gd/stepp
    dp.sir.index <- peeling(m = sqrt(rowSums(dp.sir^2)) , sig = sigma.e1, s = s, epsilon = epsilon/TT/6,delta = delta/TT/6,k=K)
    dp.sir.index
    dp.sir[dp.sir.index,] <- dp.sir[dp.sir.index,] + E1
    dp.sir[-dp.sir.index,] <- 0
  }
  error.dp.sir <- norm(proj_mat(beta1)-proj_mat(dp.sir),type = "F")
  
  return(c(error.sir,error.ini,error.ray,error.dp.sir,K))
}

# --------------------------------------------
# Linear Model Simulation
# --------------------------------------------

# Set up parallel computing

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(1004)
n <- 1000
p <- 1000
sample_cor <- ar_cov(p,0.5,0.25)
result1 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 3 ,Cn = 2, lambda = 1,stepp = 100,s=5)
}
stopCluster(cl)
# Save results to file
write.csv(result1,"result/hd_linear_n1k_p1k.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(1004)
n <- 1000
p <- 2000
result2 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 3 ,Cn = 2, lambda = 1,stepp = 110,s=5)
}
stopCluster(cl)
write.csv(result2,"result/hd_linear_n1k_p2k.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 2000
p <- 1000
result3 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 3 ,Cn = 2, lambda = 1,stepp = 100,s=5)
}
stopCluster(cl)
write.csv(result3,"result/hd_linear_n2k_p1k.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 2000
p <- 2000
result4 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 3 ,Cn = 2, lambda = 1,stepp = 110,s=5)
}
stopCluster(cl)
write.csv(result4,"result/hd_linear_n2k_p2k.csv")

# --------------------------------------------
# Exp Model Simulation
# --------------------------------------------

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(19951004)
n <- 1000
p <- 1000
result1 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 5 ,Cn = 2, lambda = 1,stepp = 100,s=5)
}
stopCluster(cl)
write.csv(result1,"result/hd_exp_n1k_p1k.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 1000
p <- 2000
result2 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 4 ,Cn = 2, lambda = 1,stepp = 100,s=5)
}
stopCluster(cl)
write.csv(result2,"result/hd_exp_n1k_p2k.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 2000
p <- 1000
result3 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 4 ,Cn = 2, lambda = 1,stepp = 100,s=5)
}
stopCluster(cl)
write.csv(result3,"result/hd_exp_n2k_p1k.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(1004)
n <- 2000
p <- 2000
result4 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sparse_sir(X,Y,H = 10,K = 1,TT1 = 1,TT2 = 4 ,Cn = 2, lambda = 1,stepp = 110,s=5)
}
stopCluster(cl)
write.csv(result4,"result/hd_exp_n2k_p2k.csv")

# --------------------------------------------
# Model M3 Simulation
# --------------------------------------------

set.seed(19951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 2000
p <- 2000
result1 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 10,K = 2,TT1 = 4,TT2 = 2 ,Cn = 0.01, lambda = 1,stepp = 40,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl)
write.csv(result1,"result/hd_t1_n2k_p2k.csv")

set.seed(19951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 2000
p <- 4000
result2 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 10,K = 2,TT1 = 4,TT2 = 2 ,Cn = 0.01, lambda = 1,stepp = 50,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl)
write.csv(result2,"result/hd_t1_n2k_p4k.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 4000
p <- 2000
result3 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 10,K = 2,TT1 = 4,TT2 = 3 ,Cn = 0.005, lambda = 1,stepp = 40,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl)
write.csv(result3,"result/hd_t1_n4k_p2k.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 4000
p <- 4000
result4 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 10,K = 2,TT1 = 4,TT2 = 3 ,Cn = 0.005, lambda = 1,stepp = 45,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl)
write.csv(result4,"result/hd_t1_n4k_p4k.csv")

# --------------------------------------------
# Model M4 Simulation
# --------------------------------------------

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 2000
p <- 2000
result21 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- (X%*%beta11)*exp(X%*%beta12+ ep)
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 6,K = 2,TT1 = 2,TT2 = 2 ,Cn = 0.01, lambda = 1,stepp = 20,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl)
write.csv(result21,"result/hd_t2_n2k_p2k.csv")

set.seed(19951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 2000
p <- 4000
result22 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- (X%*%beta11)*exp(X%*%beta12+ ep)
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 6,K = 2,TT1 = 2,TT2 = 2 ,Cn = 0.01, lambda = 1,stepp = 20,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl)
write.csv(result22,"result/hd_t2_n2k_p4k.csv")

set.seed(19951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 4000
p <- 2000
sample_cor <- ar_cov(p,0.5,0.25)
sample_cor_2 <- chol(sample_cor)
result23 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- (X%*%beta11)*exp(X%*%beta12+ ep)
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 6,K = 2,TT1 = 2,TT2 = 3 ,Cn = 0.003, lambda = 1,stepp = 18,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl) #1.7230000
write.csv(result23,"result/hd_t2_n4k_p2k.csv")

set.seed(19951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 4000
p <- 4000
result24 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    rho <- 0.5
    X <- matrix(0, nrow = n, ncol = p)
    X[,1] <- rnorm(n, mean = 0, sd = 0.5)
    for (j in 2:p) {
      X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
    }
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,-5),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,-5),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- (X%*%beta11)*exp(X%*%beta12+ ep)
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sparse_sir(X,Y,H = 6,K = 2,TT1 = 2,TT2 = 3 ,Cn = 0.003, lambda = 1,stepp = 22,s=5)
  }, error = function(e) {
    # Handle error here; could return NA, NULL, or a custom message
    cat("Error in iteration", mc, ":", conditionMessage(e), "\n")
    NA  # or any other value you prefer on error
  })
}
stopCluster(cl)
write.csv(result24,"result/hd_t2_n4k_p4k.csv")
