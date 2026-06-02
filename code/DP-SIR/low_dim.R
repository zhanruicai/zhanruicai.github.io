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

# Function: dp_sir
# Description:
#   Implements differentially private sliced inverse regression (DP-SIR)
#   with initialization, Rayleigh flow, and proposed optimization-based refinement.
# Inputs:
#   X       - covariate matrix
#   Y       - response vector
#   H       - number of slices
#   Cn      - penalty constant for BIC
#   KK      - privacy budget for iterative procedures
#   TT1     - number of iterations for DP Rayleigh flow
#   TT2     - number of iterations for proposed optimization
#   stepp   - step size for gradient update in optimization
#   lambda  - regularization parameter in optimization
# Output:
#   A vector of errors from different methods and the selected dimension K
# --------------------------------------------
dp_sir<-function(X, Y, H, Cn ,KK,TT1,TT2,stepp,lambda){
  # -------- Step 1: Slicing the response --------
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # Apply DP histogram method to Y
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 100)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  
  # Compute the slice mean matrix M and covariance matrix Sigma.X
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  
  # -------- Step 2: Initial estimation (DP) --------
  delta <- 1/n^{1.1}
  # Add symmetric noise to Sigma.X and M
  sigma.e1 <- 2*2*1^2/10*p/n*sqrt(2*log(2.5/delta))
  E1 <- rnorm(n = (p+1)*p/2, mean = 0, sd = sigma.e1)
  E1_noise <- matrix(data = 0, nrow = p, ncol = p)
  E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
  E1_noise <- E1_noise+t(E1_noise)
  diag(E1_noise)<- diag(E1_noise)/2
  sigma.e2 <- 2*6.01*1^2/10*p/n*sqrt(2*log(2.5/delta))
  E2 <- rnorm(n = (p+1)*p/2, mean = 0, sd = sigma.e2)
  E2_noise <- matrix(data = 0, nrow = p, ncol = p)
  E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E2
  E2_noise <- E2_noise+t(E2_noise)
  diag(E2_noise)<- diag(E2_noise)/2
  # Ensure PD covariance matrix
  Sigma.X.noise.PD <- as.matrix(nearPD(Sigma.X + E1_noise)$mat)
  # Eigen-decomposition with noise
  ini.sir.estimator <- geigen(A = Martix.M + E2_noise, B = Sigma.X.noise.PD, symmetric = T)
  #BIC
  order.lambda <- ini.sir.estimator$values[order(ini.sir.estimator$values,decreasing = T)]
  K = dp_bic(order.lambda = order.lambda, Cn = Cn, HH = length(unique(tildey)))
  ini.sir.number <- order(ini.sir.estimator$values,decreasing = T)[1:K]
  ini.sir.vector <- ini.sir.estimator$vectors[,ini.sir.number]
  error.ini <- norm(proj_mat(beta1,Sigma.X)-proj_mat(ini.sir.vector,Sigma.X),type = "F")
  # -------- Step 3: Vanilla SIR (non-private) --------
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.sir <- norm(proj_mat(beta1,Sigma.X)-proj_mat(sir.vector,Sigma.X),type = "F")
  # -------- Step 4: DP Rayleigh Flow Iteration --------
  epsilon <- KK
  delta <- 1/n^{1.1}
  TT <- TT1
  if (K==1) {dp.rf <- ini.sir.vector} else {dp.rf <- ini.sir.vector[,1]}
  for (T.iter in 1:TT) {
    sigma.e1 <- 2*2*1^2*p/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
    sigma.e2 <- 2*6.01*1^2*p/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
    # Noise addition
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
    # Rayleigh quotient update
    rhot <- as.numeric(t(dp.rf)%*%(Martix.M + E2_noise)%*%dp.rf/(t(dp.rf)%*%as.matrix( Sigma.X.noise.PD)%*%dp.rf))
    Ct <- diag(p)+((Martix.M + E2_noise)-rhot* Sigma.X.noise.PD)/rhot
    dp.rf <- Ct%*%dp.rf/norm(Ct%*%dp.rf,type = "2")
  }
  error.ray <- norm(proj_mat(beta1,Sigma.X)-proj_mat(dp.rf,Sigma.X),type = "F")
  
  # -------- Step 5: Proposed Optimization Approach (DP Gradient Descent) --------
  TT <- TT2
  if (K==1) {dp.sir <- ini.sir.vector*sqrt(max(ini.sir.estimator$values)/lambda+1)}
  if (K>1) {dp.sir <- ini.sir.vector%*%diag(sqrt(ini.sir.estimator$values[order(ini.sir.estimator$values,decreasing = T)[1:K]]/lambda+1))}
  #dp.sir <- sir.vector
  for (T.iter in 1:TT) {
    sigma.e1 <-  (6.01*3+lambda*2*3)*1/n*sqrt(p*K)*sqrt(2*log(1.25*TT/delta))/epsilon*TT/stepp
    gd <- -(Martix.M)%*%dp.sir+lambda*Sigma.X%*%dp.sir%*%(t(dp.sir)%*%Sigma.X%*%dp.sir-diag(K))
    E1 <- rnorm(n = p*K,mean = 0, sd = sigma.e1)
    dim(E1) <- c(p,K)
    dp.sir <- dp.sir - gd/stepp + E1
  }
  error.dp.sir <- norm(proj_mat(beta1,Sigma.X)-proj_mat(dp.sir,Sigma.X),type = "F")
  
  # -------- Return all errors and selected dimension --------
  return(c(error.sir,error.ini,error.ray,error.dp.sir,K))
}
# --------------------------------------------
# Linear Model Simulation
# --------------------------------------------

# Set up parallel computing
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 20000
p <- 15
# Generate AR(1) covariance matrix
sample_cor <- ar_cov(p,0.5,0.25)
# Run 1000 Monte Carlo simulations in parallel
result1 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  # Generate covariate matrix X from multivariate normal distribution
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  # Truncate values to [-1.5, 1.5]
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  # Create sparse coefficient vector
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  # Generate response with Gaussian noise
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  # Run DP-SIR method
  dp_sir(X,Y,H = 20,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
# Save results to file
write.csv(result1,"result/ld_linear_n2w_p15.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 20000
p <- 30
sample_cor <- ar_cov(p,0.5,0.25)
result2 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 20,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
write.csv(result2,"result/ld_linear_n2w_p30.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 40000
p <- 15
sample_cor <- ar_cov(p,0.5,0.25)
result3 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 20,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
write.csv(result3,"result/ld_linear_n4w_p15.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 40000
p <- 30
sample_cor <- ar_cov(p,0.5,0.25)
result4 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 20,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
write.csv(result4,"result/ld_linear_n4w_p30.csv")

# --------------------------------------------
# Exp Model Simulation
# --------------------------------------------
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 20000
p <- 15
sample_cor <- ar_cov(p,0.5,0.25)
result1 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
write.csv(result1,"result/ld_exp_n2w_p15.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 20000
p <- 30
sample_cor <- ar_cov(p,0.5,0.25)
result2 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
write.csv(result2,"result/ld_exp_n2w_p30.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 40000
p <- 15
sample_cor <- ar_cov(p,0.5,0.25)
result3 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
write.csv(result3,"result/ld_exp_n4w_p15.csv")

cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
set.seed(951004)
n <- 40000
p <- 30
sample_cor <- ar_cov(p,0.5,0.25)
result4 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- exp(X%*%beta1) + ep
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,KK = 1,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=2)
}
stopCluster(cl)
write.csv(result4,"result/ld_exp_n4w_p30.csv")

# --------------------------------------------
# Model M3 Simulation
# --------------------------------------------
set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 30000
p <- 10
sample_cor <- ar_cov(p,0.5,0.25)
result1 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  tryCatch({
    X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
    X[which(X>1.5)]<-1.5
    X[which(X< -1.5)]<--1.5
    beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
    beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
    ep <- rnorm(n = n, mean = 0, sd = 1)
    Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
    beta1<-cbind(beta11,beta12)
    X[which(X>1)]<-1
    X[which(X< -1)]<--1
    dp_sir(X,Y,H = 30,KK = 2,TT1 = 1,TT2=2,stepp = 1,lambda=1,Cn=0.01)}, error = function(e) {
      # Handle the error and return NA or a default value
      return(NA)
    })
}
stopCluster(cl)
write.csv(result1,"result/ld_t1_n3w_p10.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 30000
p <- 15
sample_cor <- ar_cov(p,0.5,0.25)
result2 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,KK = 2,TT1 = 1,TT2=2,stepp = 1,lambda=1,Cn=0.01)
}
stopCluster(cl)
write.csv(result2,"result/ld_t1_n3w_p15.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 50000
p <- 10
sample_cor <- ar_cov(p,0.5,0.25)
result3 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,KK = 2,TT1 = 1,TT2=2,stepp = 1,lambda=1,Cn=0.008)
}
stopCluster(cl)
write.csv(result3,"result/ld_t1_n5w_p10.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 50000
p <- 15
sample_cor <- ar_cov(p,0.5,0.25)
result4 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,KK = 2,TT1 = 1,TT2=3,stepp = 1.5,lambda=1,Cn=0.01)
}
stopCluster(cl)
write.csv(result4,"result/ld_t1_n5w_p15.csv")

# --------------------------------------------
# Model M4 Simulation
# --------------------------------------------

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 30000
p <- 10
sample_cor <- ar_cov(p,0.5,0.25)
result1 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- (X%*%beta11)*exp(X%*%beta12+ ep) 
  beta1<-cbind(beta11,beta12)
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,K = 2,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=0.01)
}
stopCluster(cl)
write.csv(result1,"result/ld_t2_n3w_p10.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 30000
p <- 15
sample_cor <- ar_cov(p,0.5,0.25)
result2 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- (X%*%beta11)*exp(X%*%beta12+ ep) 
  beta1<-cbind(beta11,beta12)
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,K = 2,TT1 = 1,TT2=1,stepp = 1,lambda=1,Cn=0.01)
}
stopCluster(cl)
write.csv(result2,"result/ld_t2_n3w_p15.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 50000
p <- 10
sample_cor <- ar_cov(p,0.5,0.25)
result3 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- (X%*%beta11)*exp(X%*%beta12+ ep) 
  beta1<-cbind(beta11,beta12)
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,K = 2,TT1 = 1,TT2=2,stepp = 1,lambda=1,Cn=0.01)
}
stopCluster(cl)
write.csv(result3,"result/ld_t2_n5w_p10.csv")

set.seed(951004)
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
n <- 50000
p <- 15
sample_cor <- ar_cov(p,0.5,0.25)
result4 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- (X%*%beta11)*exp(X%*%beta12+ ep) 
  beta1<-cbind(beta11,beta12)
  X[which(X>1)]<-1
  X[which(X< -1)]<--1
  dp_sir(X,Y,H = 30,K = 2,TT1 = 1,TT2=2,stepp = 1,lambda=1,Cn=0.01)
}
stopCluster(cl)
write.csv(result4,"result/ld_t2_n5w_p15.csv")
