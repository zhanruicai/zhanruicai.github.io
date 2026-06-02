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
library(ggplot2)
library(patchwork)

# Load all functions and variables defined in functions.R into the current R session
source("functions.R")

# --------------------------------------------
# Function: Compare DP histogram vs ground truth in SIR
# Input: X - covariate;
#        Y - response; 
#        H - slices; 
#        m - histogram bins; 
#        epsilon - privacy budget; 
#        K - dimension
# Output: A vector of [eigenvalue gap, subspace error]
# --------------------------------------------

dp_hist_vs_groud <- function(X, Y,H = 10,m=20,epsilon=1,K=1){
  n <- length(Y)
  p <- ncol(X)
  # Get noisy slicing points using DP histogram
  slicey <- dp_hist(Y = Y, H = H, epsilon = epsilon, m = m)
  # Assign slice indices to Y values
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  Sigma.X <- cov(X)
  # Compute SIR kernel matrix (second moment of slice means)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # Compute eigen decomposition of SIR matrix
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  # Compute eigenvalue gap for K-th direction
  es <- sir.estimator$values[p-(K-1)]-sir.estimator$values[p-(K-1)-1]
  # Extract leading K eigenvectors
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  # Compute projection error compared to true subspace
  error.es <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")

  c(es,error.es)
}

# --------------------------------------------
# Function: Construct autoregressive (AR(1)) covariance matrix
# Input: p - number of variables; 
#        rho - correlation coefficient;
# Output: AR(1)-structured correlation matrix
# --------------------------------------------
ar_cov<-function(p,rho,trho){
  sample_cor <- matrix(data = NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      sample_cor[i,j]<-rho^{abs(i-j)}
      #if (sample_cor[i,j]<trho) {sample_cor[i,j]<-0}
    }
  }
  return(sample_cor)
}
###########################################
# Simulation: Effect of Epsilon on DP Histogram (1-dim)
###########################################

# Set number of cores for parallel computation
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
K <- 1  # Target dimension of the sufficient subspace

# Run the simulation in parallel 5000 times
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  
  # Generate autoregressive covariance matrix
  sample_cor <- ar_cov(p,0.5,0.25)
  
  # Sample size and dimension
  n <- 1000
  p <- 10
  H <- 10 # Number of slices
  
  # Generate covariates X and truncate extreme values
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  
  # Generate true coefficient vector and response
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  Sigma.X <- cov(X) # Covariance of X
  
  # Create slicing based on true quantiles
  slicey <- quantile(Y,seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  
  # Compute SIR matrix
  
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  
  # Eigen decomposition for SIR
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  tr <- sir.estimator$values[p]-sir.estimator$values[p-1]
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
  
  # Store non-private baseline results
  ep <- c(tr,error.tr)
  
  # -----------------------------
  # Private Histogram with varying epsilon
  # -----------------------------
  for (eps in seq(0.01,0.5,by = 0.01)){
    ep <- c(ep,dp_hist_vs_groud(X = X,Y=Y,epsilon = eps,H=10,m=30,K=1))
  }
  
  ep
}
# Stop parallel cluster
stopCluster(cl)

# ------------------------------------------
# Aggregate results into two data frames
# 1. Eigen-gap
# 2. Estimation error
# ------------------------------------------

# Data frame for eigen-gap comparison
data.m1.ep.eg <- data.frame(epsilon = seq(0.01,0.5,by = 0.01),private = rowMeans(result)[seq(3,102,by = 2)], non_private = rep(rowMeans(result)[1],50))
# Data frame for subspace estimation error comparison
data.m1.ep.loss <- data.frame(epsilon = seq(0.01,0.5,by = 0.01),private = rowMeans(result)[seq(4,102,by = 2)], non_private = rep(rowMeans(result)[2],50))

#############################################################
# Simulation Study: Effect of epsilon on DP-Histogram (2-dim case)
#############################################################

# Set up parallel computing environment
cl_size = 10# Number of cores to use
cl = makeCluster(cl_size)# Create cluster
registerDoParallel(cl)# Register for foreach parallel backend

# Run simulation 5000 times in parallel
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  # ----------------------------
  # Data generation setup
  # ----------------------------
  sample_cor <- ar_cov(p,0.5,0.25)# AR(1) correlation matrix
  
  n <- 1000 # Sample size
  p <- 10 # Number of predictors
  H <- 10 # Number of slices
  
  # Generate covariate matrix X from multivariate normal
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  # Define two directions beta11 and beta12 (each sparse)
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  # Response variable: Rational model with additive noise
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)# Combine into matrix
  Sigma.X <- cov(X)# Covariance of X
  # ----------------------------
  # Standard SIR with true quantiles
  # ----------------------------
  # Compute slice boundaries (true quantiles)
  slicey <- quantile(Y,seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
  # Assign slice indices to observations
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Compute between-slice covariance matrix M
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # SIR eigen decomposition
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  tr <- sir.estimator$values[p-1]-sir.estimator$values[p-2]
  sir.number <- order(sir.estimator$values,decreasing = T)[1:2]
  sir.vector <- sir.estimator$vectors[,sir.number]
  # Subspace estimation error using Frobenius norm
  error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
  # Initialize results with non-private baseline
  ep <- c(tr,error.tr)
  # ----------------------------
  # Evaluate DP histogram at various epsilon values
  # ----------------------------
  for (eps in seq(0.01,0.5,by = 0.01)){
    ep <- c(ep,dp_hist_vs_groud(X = X,Y=Y,epsilon = eps,H=10,m=30,K=2))
  }
  ep
}
stopCluster(cl)
# ----------------------------
# Summarize and organize results
# ----------------------------
# Eigen-gap vs epsilon
data.m3.ep.eg <- data.frame(epsilon = seq(0.01,0.5,by = 0.01),private = rowMeans(result)[seq(3,102,by = 2)], non_private = rep(rowMeans(result)[1],50))
# Estimation error vs epsilon
data.m3.ep.loss <- data.frame(epsilon = seq(0.01,0.5,by = 0.01),private = rowMeans(result)[seq(4,102,by = 2)], non_private = rep(rowMeans(result)[2],50))

# Reshaping data to long format for ggplot
data_long <- data.frame(
  epsilon = rep(data.m1.ep.eg$epsilon, 2),
  value = c(data.m1.ep.eg$private, data.m1.ep.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 50)
)
pm1.eg <- ggplot(data_long, aes(x = epsilon, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Epsilon (" ~ epsilon ~ ")"),
       y = expression("Eigen-gap ("~lambda[1] - lambda[2]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
data_long <- data.frame(
  epsilon = rep(data.m3.ep.eg$epsilon, 2),
  value = c(data.m3.ep.eg$private, data.m3.ep.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 50)
)
pm3.eg <- ggplot(data_long, aes(x = epsilon, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Epsilon (" ~ epsilon ~ ")"),
       y = expression("Eigen-gap ("~lambda[2] - lambda[3]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()

# Combine plots side by side with a shared legend
combined_plot.ep.eg <- (pm1.eg + pm3.eg) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.ep.eg #Figure C.2

# Reshaping data to long format for ggplot
data_long <- data.frame(
  epsilon = rep(data.m1.ep.loss$epsilon, 2),
  value = c(data.m1.ep.loss$private, data.m1.ep.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 50)
)
pm1.loss <- ggplot(data_long, aes(x = epsilon, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Epsilon (" ~ epsilon ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
data_long <- data.frame(
  epsilon = rep(data.m3.ep.loss$epsilon, 2),
  value = c(data.m3.ep.loss$private, data.m3.ep.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 50)
)
pm3.loss <- ggplot(data_long, aes(x = epsilon, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Epsilon (" ~ epsilon ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()

# Combine plots side by side with a shared legend
combined_plot.ep.loss <- (pm1.loss + pm3.loss) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.ep.loss #Figure C.3

#################################################
# Simulation Study: Effect of histogram bin size m(1-dim)
# under fixed epsilon = 1
#################################################

# Set up parallel environment
cl_size = 10                        # Number of parallel workers
cl = makeCluster(cl_size)          # Create the cluster
registerDoParallel(cl)             # Register cluster for foreach
# Run the simulation in parallel
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  ############################
  # Data generation
  ############################
  n <- 1000 # Sample size
  p <- 10# Number of predictors
  sample_cor <- ar_cov(p,0.5,0.25)# Generate AR(1) correlation matrix
  H <- 10# Number of slices
  K <- 1# Target structural dimension
  
  # Generate covariates X from multivariate normal
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  
  # Create sparse true direction beta1 and response Y
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  Sigma.X <- cov(X)
  ############################
  # Oracle (non-private) SIR
  ############################
  # Compute true quantile slice boundaries
  slicey <- quantile(Y,seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
  # Assign slices based on true quantiles
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Compute between-slice covariance matrix M
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # SIR estimation using true quantiles
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  # Eigen-gap
  tr <- sir.estimator$values[p]-sir.estimator$values[p-1]
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  # Estimation error 
  error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
  # Store non-private results
  ep <- c(tr,error.tr)
  ############################
  # Evaluate DP-Histogram with varying m (number of bins)
  ############################
  for (mm in seq(10,150,by = 5)){
    ep <- c(ep,dp_hist_vs_groud(X = X,Y=Y,epsilon = 1,H=10,m=mm,K=1))
  }
  ep
}
stopCluster(cl)
############################
# Summarize and store results
############################
# Eigen-gap vs m
data.m1.m.eg <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(3,60,by = 2)], non_private = rep(rowMeans(result)[1],29))
# Estimation error vs m
data.m1.m.loss <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(4,60,by = 2)], non_private = rep(rowMeans(result)[2],29))

##############################################################
# Simulation Study: Effect of histogram bin number m (2-dim)
# with fixed epsilon = 1 
##############################################################
# Set up parallel processing
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  ###################################
  # Data generation
  ###################################
  sample_cor <- ar_cov(p,0.5,0.25)# Generate AR covariance matrix
  n <- 1000 # Sample size
  p <- 10# Number of covariates
  H <- 10# Number of slices
  # Simulate covariates from multivariate normal
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  # Construct two true directions (sparse beta1 with two components)
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  # Simulate nonlinear response (rational function) with noise
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)
  Sigma.X <- cov(X)# Compute sample covariance matrix
  ###################################
  # Oracle (non-private) SIR
  ###################################
  # Define true slice boundaries using quantiles of Y
  slicey <- quantile(Y,seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
  # Assign each observation to a slice
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Compute the between-slice covariance matrix M
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # Generalized eigen-decomposition for SIR
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  # Compute eigen-gap between 2nd and 3rd smallest eigenvalues
  tr <- sir.estimator$values[p-1]-sir.estimator$values[p-2]
  sir.number <- order(sir.estimator$values,decreasing = T)[1:2]
  # Extract top 2 eigenvectors (K = 2)
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
  # Initialize output vector with oracle metrics
  ep <- c(tr,error.tr)
  ###################################
  # Evaluate private histogram slicing (DP-Hist)
  # under different values of m (bin count)
  ###################################
  for (mm in seq(10,150,by = 5)){
    ep <- c(ep,dp_hist_vs_groud(X = X,Y=Y,epsilon = 1,H=10,m=mm,K=2))
  }
  ep
}
stopCluster(cl)
############################
# Summarize and store results
############################
# Eigen-gap vs m
data.m3.m.eg <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(3,60,by = 2)], non_private = rep(rowMeans(result)[1],29))
# Estimation error vs m
data.m3.m.loss <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(4,60,by = 2)], non_private = rep(rowMeans(result)[2],29))

# Reshaping data to long format for ggplot
data_long <- data.frame(
  m = rep(data.m1.m.eg$m, 2),
  value = c(data.m1.m.eg$private, data.m1.m.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm1.eg <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Eigen-gap ("~lambda[1] - lambda[2]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
mm1.eg


data_long <- data.frame(
  m = rep(data.m3.m.eg$m, 2),
  value = c(data.m3.m.eg$private, data.m3.m.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm3.eg <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Eigen-gap ("~lambda[2] - lambda[3]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
mm3.eg 

# Combine plots side by side with a shared legend
combined_plot.m.eg <- (mm1.eg + mm3.eg) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.m.eg#Figure C.4

# Reshaping data to long format for ggplot
data_long <- data.frame(
  m = rep(data.m1.m.loss$m, 2),
  value = c(data.m1.m.loss$private, data.m1.m.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm1.loss <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
mm1.loss

data_long <- data.frame(
  m = rep(data.m3.m.loss$m, 2),
  value = c(data.m3.m.loss$private, data.m3.m.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm3.loss <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()

# Combine plots side by side with a shared legend
combined_plot.m.loss <- (mm1.loss + mm3.loss) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.m.loss#Figure C.5

###############################################################
# Simulation: Impact of bin number m under epsilon = 0.1
# dimension = 1
###############################################################

cl_size = 10                                  # Number of parallel workers
cl = makeCluster(cl_size)                     # Create cluster
registerDoParallel(cl)                        # Register parallel backend
# Run 5000 independent simulations
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  n <- 1000                                   # Sample size
  p <- 10                                     # Number of covariates
  H <- 10                                     # Number of slices
  K <- 1                                      # Target dimension
  sample_cor <- ar_cov(p, 0.5, 0.25)          # Generate AR covariance matrix
  
  #Generate covariates
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  # Sparse true direction
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  Sigma.X <- cov(X)
  # Oracle slicing based on empirical quantiles
  slicey <- quantile(Y,seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Compute slice mean matrix (between-slice covariance)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # Perform SIR using generalized eigen-decomposition
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  tr <- sir.estimator$values[p]-sir.estimator$values[p-1]# Eigengap
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")# Subspace error
  
  ep <- c(tr,error.tr)
  # Run private histogram slicing with various bin sizes m
  for (mm in seq(10,150,by = 5)){
    ep <- c(ep,dp_hist_vs_groud(X = X,Y=Y,epsilon = 0.1,H=10,m=mm,K=1))
  }
  ep
}
stopCluster(cl)
# Compile result for eigen-gap and loss
data.m1.m.eg <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(3,60,by = 2)], non_private = rep(rowMeans(result)[1],29))
data.m1.m.loss <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(4,60,by = 2)], non_private = rep(rowMeans(result)[2],29))

###############################################################
# Simulation: Impact of bin number m under epsilon = 0.1
# dimension = 2
###############################################################

cl_size = 10                                  # Number of parallel workers
cl = makeCluster(cl_size)                     # Create cluster
registerDoParallel(cl)                        # Register parallel backend

result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  n <- 1000                                   # Sample size
  p <- 10                                     # Number of covariates
  H <- 10                                     # Number of slices
  K <- 2                                      # Target dimension
  sample_cor <- ar_cov(p, 0.5, 0.25)          # Generate AR covariance matrix
  #Generate Covariates
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  # Sparse true direction
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)
  Sigma.X <- cov(X)
  # Oracle slicing based on empirical quantiles
  slicey <- quantile(Y,seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Compute slice mean matrix (between-slice covariance)
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # Perform SIR using generalized eigen-decomposition
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  tr <- sir.estimator$values[p-1]-sir.estimator$values[p-2]# Eigengap
  sir.number <- order(sir.estimator$values,decreasing = T)[1:2]
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")# Subspace error
  #initial value for non-private
  ep <- c(tr,error.tr)
  # Run private histogram slicing with various bin sizes m
  for (mm in seq(10,150,by = 5)){
    ep <- c(ep,dp_hist_vs_groud(X = X,Y=Y,epsilon = 0.1,H=10,m=mm,K=2))
  }
  ep
}
stopCluster(cl)
# Compile result for eigen-gap and loss
data.m3.m.eg <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(3,60,by = 2)], non_private = rep(rowMeans(result)[1],29))
data.m3.m.loss <- data.frame(m = seq(10,150,by = 5),private = rowMeans(result)[seq(4,60,by = 2)], non_private = rep(rowMeans(result)[2],29))

# Reshaping data to long format for ggplot
data_long <- data.frame(
  m = rep(data.m1.m.eg$m, 2),
  value = c(data.m1.m.eg$private, data.m1.m.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm1.eg <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Eigen-gap ("~lambda[1] - lambda[2]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
mm1.eg

# Reshaping data to long format for ggplot
data_long <- data.frame(
  m = rep(data.m3.m.eg$m, 2),
  value = c(data.m3.m.eg$private, data.m3.m.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm3.eg <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Eigen-gap ("~lambda[2] - lambda[3]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
mm3.eg 

# Combine plots side by side with a shared legend
combined_plot.m.eg <- (mm1.eg + mm3.eg) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.m.eg#Figure C.6

# Reshaping data to long format for ggplot
data_long <- data.frame(
  m = rep(data.m1.m.loss$m, 2),
  value = c(data.m1.m.loss$private, data.m1.m.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm1.loss <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
mm1.loss
# Reshaping data to long format for ggplot
data_long <- data.frame(
  m = rep(data.m3.m.loss$m, 2),
  value = c(data.m3.m.loss$private, data.m3.m.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 29)
)
mm3.loss <- ggplot(data_long, aes(x = m, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Number of bins (" ~ m ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()

# Combine plots side by side with a shared legend
combined_plot.m.loss <- (mm1.loss + mm3.loss) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.m.loss#Figure C.7

# ------------------------------------------------------------------
# Function: dp_hist_vs_groud_both
# Purpose:
#   Compare differentially private (DP) histogram vs. true quantile histogram
#   in terms of eigengap and subspace estimation error from SIR.
#
# Inputs:
#   - X: n x p covariate matrix
#   - Y: n-dimensional response vector
#   - H: number of slices
#   - m: number of histogram bins
#   - epsilon: DP privacy budget
#   - K: target structural dimension
#
# Output:
#   - A vector of [eigenvalue gap, subspace error]
# ------------------------------------------------------------------

dp_hist_vs_groud_both <- function(X, Y,H = 10,m=20,epsilon=1,K=1){
  n <- length(Y)
  p <- ncol(X)
  
  #####################
  # --- DP Histogram ---#
  #####################
  
  # Step 1: Get differentially private slice boundaries
  slicey <- dp_hist(Y = Y, H = H, epsilon = epsilon, m = m)
  # Step 2: Assign observations to slices based on DP quantiles
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Step 3: Compute covariance matrix of X
  Sigma.X <- cov(X)
  # Step 4: Construct between-slice covariance matrix M
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # Step 5: Perform SIR and extract eigengap and subspace error
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  es <- sir.estimator$values[p-(K-1)]-sir.estimator$values[p-(K-1)-1]
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.es <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
  ###########################
  # --- Non-private Histogram ---#
  ###########################
  
  Sigma.X <- cov(X)
  # Step 6: Use true quantiles for slicing
  slicey <- quantile(Y,seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
  # Step 7: Re-assign slices using empirical quantiles
  tildey <- rep(1,n)
  for (i in 1:(H-1)) {
    tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
  }
  # Step 8: Reconstruct between-slice covariance matrix M
  Martix.M <- matrix(data = 0, nrow = p, ncol = p)
  overmean <- colMeans(X)
  for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
    if (length(which(tildey==i))==1) {slicemean <- X[which(tildey==i),]} else {
      slicemean <- colMeans(X[which(tildey==i),])}
    Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
  }
  # Step 9: Perform SIR with true slicing
  sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
  tr <- sir.estimator$values[p-(K-1)]-sir.estimator$values[p-(K-1)-1]
  sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
  sir.vector <- sir.estimator$vectors[,sir.number]
  error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
  ###########################
  # Output: Compare both methods
  ###########################
  c(tr,error.tr,es,error.es)
}
###########################################
# Experiment: Effect of number of slices H under fixed epsilon(dim = 1)
# Goal:
#   Evaluate how the number of slices (H) impacts the eigengap and subspace estimation error
#   under both differentially private and non-private slicing.
# Parameters:
#   - Fixed epsilon = 0.1
#   - Vary H from 2 to 20
#   - Fixed number of histogram bins (m = 50)
###########################################

# Set up parallel computation with 10 cores
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)

# Run simulations in parallel over 5000 Monte Carlo repetitions
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  sample_cor <- ar_cov(p,0.5,0.25)# AR covariance structure
  # Data generation setup
  p <- 10           # Number of covariates
  n <- 1000         # Sample size
  K <- 1            # Structural dimension
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  # Response variable
  Y <- X%*%beta1 + ep
  # Loop over number of slices H (from 2 to 20)
  ep <- c()
  for (HH in seq(2,20,by = 1)){
    ep <- c(ep,dp_hist_vs_groud_both(X = X,Y=Y,epsilon = 0.1,H=HH,m=50,K=1))
  }
  ep
}
stopCluster(cl)
# Organize results into two data frames for plotting or reporting
data.m1.H.eg <- data.frame(H = seq(2,20,by = 1),private = rowMeans(result)[seq(3,76,by = 4)], non_private = rowMeans(result)[seq(1,76,by = 4)])
data.m1.H.loss <- data.frame(H = seq(2,20,by = 1),private = rowMeans(result)[seq(4,76,by = 4)], non_private = rowMeans(result)[seq(2,76,by = 4)])
###########################################
# Experiment: Effect of number of slices H under fixed epsilon(dim = 2)
# Goal:
#   Evaluate how the number of slices (H) impacts the eigengap and subspace estimation error
#   under both differentially private and non-private slicing.
# Parameters:
#   - Fixed epsilon = 0.1
#   - Vary H from 2 to 20
#   - Fixed number of histogram bins (m = 50)
###########################################

# Set up parallel computation with 10 cores
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)

result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  sample_cor <- ar_cov(p,0.5,0.25)
  # Parameters
  p <- 10       # Number of features
  n <- 1000     # Sample size
  H <- 10       # Number of slices (initial, varies below)
  K <- 2        # True structural dimension
  # Generate sample covariates
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  # Generate two sparse direction vectors
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  # Generate nonlinear response using t1-type model structure
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)
  
  # Loop over H = 2 to 20 and evaluate both private and non-private performance
  ep <- c()
  for (HH in seq(2,20,by = 1)){
    ep <- c(ep,dp_hist_vs_groud_both(X = X,Y=Y,epsilon = 0.1,H=HH,m=50,K=2))
  }
  ep
}
stopCluster(cl)
# Organize result into data frames for visualization/comparison
data.m3.H.eg <- data.frame(H = seq(2,20,by = 1),private = rowMeans(result)[seq(3,76,by = 4)], non_private = rowMeans(result)[seq(1,76,by = 4)])
data.m3.H.loss <- data.frame(H = seq(2,20,by = 1),private = rowMeans(result)[seq(4,76,by = 4)], non_private = rowMeans(result)[seq(2,76,by = 4)])

# Reshaping data to long format for ggplot
data_long <- data.frame(
  H = rep(data.m1.H.eg$H, 2),
  value = c(data.m1.H.eg$private, data.m1.H.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 19)
)
Hm1.eg <- ggplot(data_long, aes(x = H, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = "H",
       y = expression("Eigen-gap ("~lambda[1] - lambda[2]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
Hm1.eg
# Reshaping data to long format for ggplot
data_long <- data.frame(
  H = rep(data.m3.H.eg$H, 2),
  value = c(data.m3.H.eg$private, data.m3.H.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 19)
)
Hm3.eg <- ggplot(data_long, aes(x = H, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = "H",
       y = expression("Eigen-gap ("~lambda[2] - lambda[3]~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
Hm3.eg
# Combine plots side by side with a shared legend
combined_plot.H.eg <- (Hm1.eg + Hm3.eg) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.H.eg#Figure C.8

# Reshaping data to long format for ggplot
data_long <- data.frame(
  H = rep(data.m1.H.loss$H, 2),
  value = c(data.m1.H.loss$private, data.m1.H.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 19)
)
Hm1.loss <- ggplot(data_long, aes(x = H, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = "H",
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
Hm1.loss
# Reshaping data to long format for ggplot
data_long <- data.frame(
  H = rep(data.m3.H.loss$H, 2),
  value = c(data.m3.H.loss$private, data.m3.H.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 19)
)
Hm3.loss <- ggplot(data_long, aes(x = H, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = "H",
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
Hm3.loss

# Combine plots side by side with a shared legend
combined_plot.H.loss <- (Hm1.loss + Hm3.loss) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.H.loss#Figure C.9
###############################################################
# Experiment: Sample Size Impact on DP vs Non-DP Performance(dim=1)
# Model: Single-index linear model
# Goal:
#   - Investigate how increasing Epsilon:
#       1. Eigengap (separation between relevant and irrelevant eigenspaces)
#       2. Subspace estimation error
#   - Compare the performance between private and non-private SIR
#   - Use adaptive histogram bin size (m ~ n^{1/3})
###############################################################

# Set up parallel backend with 10 workers
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  # Global parameters
  n <- 4000      # Total sample size
  p <- 10        # Number of covariates
  K <- 1         # True structural dimension
  H <- 10        # Number of slices
  # Generate covariance matrix and covariates
  sample_cor <- ar_cov(p, 0.5, 0.25)
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)]<-1.5
  X[which(X< -1.5)]<--1.5
  # Generate sparse direction and response
  beta1 <- c(runif(n=2,-10,10),rep(0,p-2))
  dim(beta1)<-c(p,1)
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- X%*%beta1 + ep
  
  ep <- c()
  # Loop through increasing eps
  for (eps in seq(500,4000,by = 50)){
    ## ---------- Non-private SIR Estimation ---------- ##
    # Compute slice boundaries based on empirical quantiles
    slicey <- quantile(Y[1:eps],seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
    tildey <- rep(1,eps)
    for (i in 1:(H-1)) {
      tildey[which(Y[1:eps]<=slicey[i])] <- tildey[which(Y[1:eps]<=slicey[i])] + 1
    }
    # Compute SIR matrix and covariance matrix using first eps samples
    Martix.M <- matrix(data = 0, nrow = p, ncol = p)
    overmean <- colMeans(X[1:eps,])
    Sigma.X <- cov(X[1:eps,])
    for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
      if (length(which(tildey==i))==1) {slicemean <- (X[1:eps,])[which(tildey==i),]} else {
        slicemean <- colMeans((X[1:eps,])[which(tildey==i),])}
      Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/eps)
    }
    # Compute eigen-decomposition
    sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
    tr <- sir.estimator$values[p]-sir.estimator$values[p-1]
    sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
    sir.vector <- sir.estimator$vectors[,sir.number]
    error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
    ## ---------- Private SIR Estimation ---------- ##
    # Histogram bin size m scales with n^{1/3}
    ep <- c(ep,dp_hist_vs_groud(X = X[1:eps,],Y=Y[1:eps],epsilon = 0.1,H=10,m=floor(8*eps^{1/3}),K=1),sir.estimator$values[p]-sir.estimator$values[p-1],error.tr)
  }
  ep
}
stopCluster(cl)
## Post-process results into two data.frames for plotting / analysis
data.m1.n.eg <- data.frame(n = seq(500,4000,by = 50),private = rowMeans(result)[seq(1,284,by = 4)], non_private = rowMeans(result)[seq(3,284,by = 4)])
data.m1.n.loss <- data.frame(n = seq(500,4000,by = 50),private = rowMeans(result)[seq(2,284,by = 4)], non_private = rowMeans(result)[seq(4,284,by = 4)])
###############################################################
# Experiment: Sample Size Impact on DP vs Non-DP Performance(dim=2)
# Model: Single-index linear model
# Goal:
#   - Investigate how increasing Epsilon:
#       1. Eigengap (separation between relevant and irrelevant eigenspaces)
#       2. Subspace estimation error
#   - Compare the performance between private and non-private SIR
#   - Use adaptive histogram bin size (m ~ n^{1/3})
###############################################################

# Set up parallel computing environment with 10 workers
cl_size = 10
cl = makeCluster(cl_size)
registerDoParallel(cl)
result <-  foreach(mc = 1:5000, .combine = cbind,.packages = c("geigen","VGAM","Matrix","MASS")) %dopar% {
  # Simulation parameters
  n <- 5000           # Maximum sample size
  p <- 10             # Number of predictors
  H <- 9              # Number of slices
  K <- 2              # Structural dimension
  # Generate autoregressive covariance and covariates
  sample_cor <- ar_cov(p, 0.5, 0.25)
  X <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_cor*0.25)
  X[which(X>1.5)] <- 1.5
  X[which(X< -1.5)]<- -1.5
  # Generate response
  beta11 <- c(runif(n=2,-10,10),rep(0,p-2))
  beta12 <- c(runif(n=2,-10,10),rep(0,p-2))
  ep <- rnorm(n = n, mean = 0, sd = 1)
  Y <- 25*(X%*%beta11)/(1+(X%*%beta12+1)^2) + 0.1*ep
  beta1<-cbind(beta11,beta12)
  ep <- c()
  # Loop over growing epsilon
  for (eps in seq(500,5000,by = 50)){
    ########## Non-Private SIR (True Quantile Slicing) ##########
    
    # Step 1: Compute slicing based on true quantiles
    slicey <- quantile(Y[1:eps],seq(from = 0, to = 1,length.out=H+1)[-c(1,H+1)])
    tildey <- rep(1,eps)
    for (i in 1:(H-1)) {
      tildey[which(Y[1:eps]<=slicey[i])] <- tildey[which(Y[1:eps]<=slicey[i])] + 1
    }
    # Step 2: Compute between-slice covariance matrix
    Martix.M <- matrix(data = 0, nrow = p, ncol = p)
    overmean <- colMeans(X[1:eps,])
    Sigma.X <- cov(X[1:eps,])
    for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
      if (length(which(tildey==i))==1) {slicemean <- (X[1:eps,])[which(tildey==i),]} else {
        slicemean <- colMeans((X[1:eps,])[which(tildey==i),])}
      Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/eps)
    }
    # Step 3: SIR eigen decomposition
    sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
    sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
    sir.vector <- sir.estimator$vectors[,sir.number]
    error.tr <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
    tr <- sir.estimator$values[10-1]-sir.estimator$values[10-2]
    ########## Private SIR (DP Histogram Slicing) ##########
    
    # Step 1: Estimate DP quantile boundaries using histogram
    slicey <- dp_hist(Y = Y[1:eps], H = 9, epsilon = 0.1, m = floor(5.3*eps^{1/3}))
    tildey <- rep(1,eps)
    for (i in 1:(H-1)) {
      tildey[which(Y[1:eps]<=slicey[i])] <- tildey[which(Y[1:eps]<=slicey[i])] + 1
    }
    # Step 2: Recompute SIR matrix using DP-based slicing
    Martix.M <- matrix(data = 0, nrow = p, ncol = p)
    overmean <- colMeans(X[1:eps,])
    for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
      if (length(which(tildey==i))==1) {slicemean <- (X[1:eps,])[which(tildey==i),]} else {
        slicemean <- colMeans((X[1:eps,])[which(tildey==i),])}
      Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/eps)
    }
    # Step 3: Eigen decomposition with DP slicing
    sir.estimator <- geigen(A = Martix.M, B = Sigma.X, symmetric = T)
    es <- sir.estimator$values[p-(K-1)]-sir.estimator$values[p-(K-1)-1]
    sir.number <- order(sir.estimator$values,decreasing = T)[1:K]
    sir.vector <- sir.estimator$vectors[,sir.number]
    error.es <- norm(proj_mat(beta1)-proj_mat(sir.vector),type = "F")
    
    # Store all metrics: (DP eigengap, DP error, non-DP eigengap, non-DP error)
    ep <- c(ep,es,error.es,tr,error.tr)
  }
  ep
}
stopCluster(cl)
# Create data frames for plotting/analysis
data.m3.n.eg <- data.frame(n = seq(500,5000,by = 50),private = rowMeans(result)[seq(1,364,by = 4)], non_private = rowMeans(result)[seq(3,364,by = 4)])
data.m3.n.loss <- data.frame(n = seq(500,5000,by = 50),private = rowMeans(result)[seq(2,364,by = 4)], non_private = rowMeans(result)[seq(4,364,by = 4)])

# Reshaping data to long format for ggplot
data_long <- data.frame(
  n = rep(data.m1.n.eg$n, 2),
  value = c(data.m1.n.eg$private, data.m1.n.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 71)
)
pm1.n.eg <- ggplot(data_long, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Sample size (" ~ n ~ ")"),
       y = expression("Eigen-gap ("~lambda[1] - lambda[2]~ ")"),
       color = "Group",
       linetype = "Group") + ylim(c(0.65,0.85))+
  theme_minimal()
pm1.n.eg

# Reshaping data to long format for ggplot
data_long <- data.frame(
  n = rep(data.m1.n.loss$n, 2),
  value = c(data.m1.n.loss$private, data.m1.n.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 71)
)
pm1.n.loss <- ggplot(data_long, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M1):k=1",
       x = expression("Sample size (" ~ n ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
pm1.n.loss
# Reshaping data to long format for ggplot
data_long <- data.frame(
  n = rep(data.m3.n.eg$n, 2),
  value = c(data.m3.n.eg$private, data.m3.n.eg$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 91)
)
pm3.n.eg <- ggplot(data_long, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Sample size (" ~ n ~ ")"),
       y = expression("Eigen-gap ("~lambda[2] - lambda[3]~ ")"),
       color = "Group",
       linetype = "Group") + ylim(c(0.075,0.221)) + 
  theme_minimal()
pm3.n.eg
# Reshaping data to long format for ggplot
data_long <- data.frame(
  n = rep(data.m3.n.loss$n, 2),
  value = c(data.m3.n.loss$private, data.m3.n.loss$non_private),
  group = rep(c("DP Histogram", "Non-Private Quantile"), each = 91)
)
pm3.n.loss <- ggplot(data_long, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = c("DP Histogram" = "black", "Non-Private Quantile" = "red")) +
  scale_linetype_manual(values = c("DP Histogram" = "solid", "Non-Private Quantile" = "dashed")) +
  labs(title = "(M3):k=2",
       x = expression("Sample size (" ~ n ~ ")"),
       y = expression("Loss ("~L~ ")"),
       color = "Group",
       linetype = "Group") +
  theme_minimal()
pm3.n.loss

combined_plot.n.eg <- (pm1.n.eg + pm3.n.eg) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.n.eg#Figure C.10

combined_plot.n.loss <- (pm1.n.loss + pm3.n.loss) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot.n.loss#Figure C.11
