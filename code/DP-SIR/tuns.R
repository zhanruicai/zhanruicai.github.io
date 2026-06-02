# --- Load required libraries ---
library(geigen)
library(VGAM)
library(Matrix)
library(foreach)
library(parallel)
library(doParallel)
library(MASS)

# Load all functions and variables defined in functions.R into the current R session
source("functions.R")

# --- Main Function: Differentially Private Sparse SIR Estimation ---
dp_sparse_sir_simple<-function(X, Y, H, K, s,TT2,stepp,lambda){
  n <- dim(X)[1]
  p <- dim(X)[2]
  # Discretize response Y using private histogram
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Construct SIR covariance matrix
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  
  # --- Initial Estimation via Noisy Generalized Eigen Decomposition ---
  epsilon <- 2
  delta <- 1/n^{1.1}
  # Private selection of pindex (important variables)
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = s,epsilon = 3.5*epsilon, delta = delta,k=K)
  # Add noise to numerator and denominator matrices
  sigma.e1 <- 2*s/n*sqrt(2*log(2.5/delta))/epsilon
  E1 <- rnorm(n = (s+1)*s/2, mean = 0, sd = sigma.e1)
  E1_noise <- matrix(data = 0, nrow = s, ncol = s)
  E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
  E1_noise <- E1_noise+t(E1_noise)
  diag(E1_noise)<- diag(E1_noise)/2
  sigma.e2 <- 2*s/n*sqrt(2*log(2.5/delta))/epsilon
  E2 <- rnorm(n = (s+1)*s/2, mean = 0, sd = sigma.e2)
  E2_noise <- matrix(data = 0, nrow = s, ncol = s)
  E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E2
  E2_noise <- E2_noise+t(E2_noise)
  diag(E2_noise)<- diag(E2_noise)/2
  # Compute private generalized eigenvalue decomposition
  Sigma.X.noise.PD <- as.matrix(nearPD(Sigma.X[pindex,pindex] + E1_noise)$mat)
  ini.sir.estimator <- geigen(A = Martix.M[pindex,pindex] + E2_noise, B = Sigma.X.noise.PD, symmetric = T)
  ini.sir.number <- order(ini.sir.estimator$values,decreasing = T)[1:K]
  ini.sir.vector <- ini.sir.estimator$vectors[,ini.sir.number]
  ini.sir.vector.complete <-matrix(data = 0, nrow = p, ncol = K)
  ini.sir.vector.complete[pindex,] <- ini.sir.vector
  # Initial estimation error (compared to true direction beta1)
  error.ini <- norm(proj_mat(beta1)-proj_mat((ini.sir.vector.complete)),type = "F")

  # --- Proposed Method: Penalized DP-SIR ---
  if (K==1) {dp.sir <- ini.sir.vector.complete*sqrt(max(ini.sir.estimator$values)/lambda+1)}
  if (K>1) {dp.sir <- ini.sir.vector.complete%*%diag(sqrt(ini.sir.estimator$values[order(ini.sir.estimator$values,decreasing = T)[1:K]]/lambda+1))}
  dim(dp.sir)<-c(p,K)
  TT <- TT2# Number of iterations
  for (T.iter in 1:TT) {
    sigma.e1 <- (6.01*4+lambda*2*4)/n*sqrt(K)/stepp
    gd <- -Martix.M%*%dp.sir+lambda*Sigma.X%*%dp.sir%*%(t(dp.sir)%*%Sigma.X%*%dp.sir-diag(K))
    gd <- gd
    # Add noise and perform private gradient update with sparsity
    E1 <- rnorm(n = s*K, mean = 0, sd = 6/5*sigma.e1*sqrt(s)*sqrt(2*log(1.25*TT*6/5/delta))/epsilon*TT)
    dp.sir <- dp.sir - gd/stepp
    dp.sir.index <- peeling(m = sqrt(rowSums(dp.sir^2)) , sig = sigma.e1, s = s, epsilon = epsilon/TT/6,delta = delta/TT/6,k=K)
    dp.sir.index
    dp.sir[dp.sir.index,] <- dp.sir[dp.sir.index,] + E1
    dp.sir[-dp.sir.index,] <- 0
  }
  # Final estimation error
  error.dp.sir <- norm(proj_mat(beta1)-proj_mat(dp.sir),type = "F")
  
  return(c(error.ini,error.dp.sir))
}
# --- Differentially Private Cross-Validation for SIR Subspace Selection ---
dp_cv<-function(Xtr, Ytr,Xte,Yte,H, K, srange){
  ntr <- dim(Xtr)[1]
  ptr <- dim(Xtr)[2]
  # --- Compute covariance and between-slice covariance matrix on training data ---
  Sigma.X <- cov(Xtr)
  Martix.M <- matrix(data = 0, nrow = ptr, ncol = ptr)
  overmean <- colMeans(Xtr)
  for (i in unique(Ytr)[order(unique(Ytr),decreasing = F)]) {
       if (length(which(Ytr==i))==1) {slicemean <- Xtr[which(Ytr==i),]} else {
            slicemean <- colMeans(Xtr[which(Ytr==i),])}
         Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(Ytr==i))/ntr)
       }
     
  ##### --- Add Differential Privacy Noise to the Training Matrices ---
    sigma.e1 <- ptr/n*sqrt(2*log(2.5/delta))/epsilon
    E1 <- rnorm(n = (ptr+1)*ptr/2, mean = 0, sd = sigma.e1)
    E1_noise <- matrix(data = 0, nrow = ptr, ncol = ptr)
    E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
    E1_noise <- E1_noise+t(E1_noise)
    diag(E1_noise)<- diag(E1_noise)/2
    sigma.e2 <- ptr/n*sqrt(2*log(2.5/delta))/epsilon
    E2 <- rnorm(n = (ptr+1)*ptr/2, mean = 0, sd = sigma.e2)
    E2_noise <- matrix(data = 0, nrow = ptr, ncol = ptr)
    E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E1
    E2_noise <- E2_noise+t(E2_noise)
    diag(E2_noise)<- diag(E2_noise)/2
    # Ensure the perturbed covariance matrix is positive definite
    Sigma.X.noise.PD <- as.matrix(nearPD(Sigma.X + E2_noise)$mat)
    
    # --- Compute test-set between-slice covariance matrix ---
    t.M <- matrix(data = 0, nrow = ptr, ncol = ptr)
    overmean <- colMeans(Xte)
    for (i in unique(Yte)[order(unique(Yte),decreasing = F)]) {
      if (length(which(Yte==i))==1) {slicemean <- Xte[which(Yte==i),]} else {
        slicemean <- colMeans(Xte[which(Yte==i),])}
        t.M  <- t.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(Yte==i))/length(Yte))
    }
    # --- For each candidate sparsity level, evaluate test-set explained variance ---
    retu <- c()
    for (s in srange) {
      ini.sir.estimator <- geigen(A = Martix.M[1:s,1:s] + E1_noise[1:s,1:s], B = Sigma.X.noise.PD[1:s,1:s], symmetric = T)
      ini.sir.number <- order(ini.sir.estimator$values,decreasing = T)[1:K]
      ini.sir.vector <- ini.sir.estimator$vectors[,ini.sir.number]
      # Compute projected test-sample variance explained by the estimated direction
      retu <- c(retu,as.numeric(ini.sir.vector%*%t.M[1:s,1:s]%*%ini.sir.vector))
    }
  return(retu)
}

# --- Setup parallel backend ---
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
# --- Simulation parameters ---
set.seed(951004)
n <- 1000
p <- 1000
# --- Begin parallel simulation (1000 Monte Carlo replications) ---
result1000 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  # --- Generate correlated covariates X (autoregressive structure) ---
  B <- 1
  rho <- 0.5
  X <- matrix(0, nrow = n, ncol = p)
  X[,1] <- rnorm(n, mean = 0, sd = 0.5)
  for (j in 2:p) {
    X[,j] <- rho * X[,j-1] + rnorm(n, mean = 0, sd = 0.5)*sqrt(1-rho^2)
  }
  # Truncate extreme values to ensure bounded input
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,-5),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  # --- Discretize Y using differentially private histogram (for slicing) ---
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # --- Compute between-slice covariance matrix (SIR) ---
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # --- Initial private feature selection using peeling ---
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  # --- 5-fold private cross-validation for sparsity tuning ---
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  # --- Exponential mechanism to select s based on CV loss ---
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  # --- Compare proposed DP-SIR using selected s and fixed s = 2 ---
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1000)
set.seed(951004)
n <- 1100
result1100 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1100)
set.seed(951004)
n <- 1200
result1200 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1200)
set.seed(951004)
n <- 1300
result1300 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1300)
set.seed(951004)
n <- 1400
result1400 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1400)
set.seed(951004)
n <- 1500
result1500 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1500)
set.seed(951004)
n <- 1600
result1600 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1600)
set.seed(951004)
n <- 1700
result1700 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1700)
set.seed(951004)
n <- 1800
result1800 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1800)
set.seed(951004)
n <- 1900
result1900 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result1900)
set.seed(951004)
n <- 2000
result2000 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2000)
set.seed(951004)
n <- 2100
p <- 2100
result2100 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2100)
set.seed(951004)
n <- 2200
result2200 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2200)
set.seed(951004)
n <- 2300
result2300 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2300)
set.seed(951004)
n <- 2400
result2400 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2400)
set.seed(951004)
n <- 2500
result2500 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2500)
set.seed(951004)
n <- 2600
result2600 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2600)
set.seed(951004)
n <- 2700
result2700 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2700)
set.seed(951004)
n <- 2800
result2800 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2800)
set.seed(951004)
n <- 2900
result2900 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result2900)
set.seed(951004)
n <- 3000
result3000 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result3000)
set.seed(951004)
n <- 3100
result3100 <-  foreach(mc = 1:1000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  #generation
  B <- 1
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
  X[which(X>1)]<- 1
  X[which(X< -1)]<- -1
  #CV
  H = 10
  smax = 8
  slicey <- dp_hist(Y = Y, H = H, epsilon = 0.1, m = 50)
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  #####initial#####
  K <- 1
  epsilon <- 1
  delta <- 1/n^{1.1}
  pindex <- peeling(m = diag(Martix.M), sig = 6.01/n, s = smax,epsilon = 3.5*epsilon, delta = delta,k=K)
  
  cvid <-  rep(1:5,each = n/5)
  cvloss <- dp_cv(Xtr = X[which(cvid!=1),pindex],Ytr = tildey[which(cvid!=1)],Xte = X[which(cvid==1),pindex],Yte = tildey[which(cvid==1)],H = 10,K = 1,srange = c(2:smax))
  
  potential <- cvloss
  pp <- exp(1*potential/(2*1.6^2*6.1/(n/5)*K))
  pp <- cumsum(pp/sum(pp))
  shat <- min(which(runif(n = 1)<=pp))+1
  c(dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=shat),dp_sparse_sir_simple(X,Y,H = 6,K = 1,TT2 = 2 ,lambda = 1,stepp = 20,s=2))
}
rowMeans(result3100)

# --- Aggregate simulation results ---
# Each resultXX00 object contains simulation results for a different sample size (from n = 1000 to n = 3000)
plotdataM1 <- cbind(rowMeans(result1000),rowMeans(result1100),rowMeans(result1200),rowMeans(result1300),rowMeans(result1400),rowMeans(result1500),rowMeans(result1600),rowMeans(result1700),rowMeans(result1800),rowMeans(result1900),rowMeans(result2000),
                    rowMeans(result2100),rowMeans(result2200),rowMeans(result2300),rowMeans(result2400),rowMeans(result2500),rowMeans(result2600),rowMeans(result2700),rowMeans(result2800),rowMeans(result2900),rowMeans(result3000))
# --- Construct tidy dataframe for ggplot ---
# plotdataM1[2, ]: DP-SIR using cross-validated s
# plotdataM1[4, ]: DP-SIR using oracle-selected s (lower bound on loss)
data_long <- data.frame(
  n = rep(seq(1,3,by=0.1), 2),
  value = c(plotdataM1[2,], plotdataM1[4,]),
  group = rep(c("Validation", "Oracle"), each = 21)
)
# --- Create performance curve: validation vs oracle ---


CVM1 <- ggplot(data_long, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("Validation" = "black", "Oracle" = "red")) +
  scale_linetype_manual(values = c("Validation" = "solid", "Oracle" = "dashed")) +
  labs(x = "Sample size/1000",
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
CVM1


