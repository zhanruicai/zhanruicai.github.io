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
library(expm)
# --------------------------------------------
# Function: Compute projection matrix
# Input:
#   x - a numeric matrix
# Output:
#   A projection matrix of the same dimension
# --------------------------------------------
proj_mat <- function(x){
  return(x%*%solve(t(x)%*%x)%*%t(x))
}
# --------------------------------------------
# Function: Estimate quantiles based on histogram CDF
# Input:
#   p      - the target quantile level (between 0 and 1)
#   cdf    - a cumulative distribution vector (same length as number of bins)
#   breaks - the bin breakpoints from a histogram
# Output:
#   The bin edge where the CDF first exceeds or equals p
# --------------------------------------------
find_quantile <- function(p, cdf, breaks) {
  index <- which(cdf >= p)[1]  # Find the first bin where CDF is greater than or equal to p
  return(breaks[index])  # Return the corresponding bin edge
}
# --------------------------------------------
# Function: Differentially Private Histogram Estimation
# Input:
#   Y       - numeric vector
#   H       - number of slices
#   epsilon - privacy budget
#   m       - number of bins for histogram
# Output:
#   A numeric vector of estimated quantile cutoffs in the original scale
# --------------------------------------------
dp_hist <- function(Y, H, epsilon, m = 0){
  if (m==0) {m=floor(length(Y)^{1/3})} # Default bin number if not specified
  TY <- atan(Y)/(pi/2) # Transformation to [-1,1] for histogram binning
  nphist <- hist(TY,breaks = seq(-1, 1, length.out = m + 1),plot = FALSE)
  npcounts <- nphist$counts
  counts <- pmax(0,npcounts + rlaplace(n = m, location = 0, scale = 2/epsilon)) # Add Laplace noise to ensure DP
  breaks <- nphist$breaks # Calculate the CDF from the histogram counts
  cdf <- cumsum(counts) / sum(counts) # Define quantiles
  quantiles <- seq(0, 1, length.out = H + 1)[c(-1,-(H+1))] # Get quantile values
  quantile_values <- sapply(quantiles, find_quantile, cdf = cdf, breaks = breaks)
  return(tan(quantile_values*(pi/2))) # Inverse transform
}
# --------------------------------------------
# Function: BIC-style criterion to select dimension K
# Inputs:
#   order.lambda - sorted eigenvalues
#   HH           - number of slices
#   Cn           - penalty term constant
# Output:
#   Optimal number of components (K)
# --------------------------------------------
dp_bic <- function(order.lambda,HH,Cn){
  G_function <- c()
  for (l in 1:HH) {
    # Cumulative eigenvalue proportion penalized by complexity term
    G_function[l] <- sum((order.lambda[1:l])^2)/sum((order.lambda)^2)-Cn*(l+1)*l/2
  }
  return(which.max(G_function))
}
# --------------------------------------------
# Function to generate an autoregressive covariance matrix
# Inputs:
#   p    - the number of variables (dimension of the matrix)
#   rho  - correlation coefficient for AR structure
# Output: 
#   AR covariance matrix
# --------------------------------------------
ar_cov<-function(p,rho){
  sample_cor <- matrix(data = NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      sample_cor[i,j]<-rho^{abs(i-j)}
      #if (sample_cor[i,j]<trho) {sample_cor[i,j]<-0}
    }
  }
  return(sample_cor)
}
# --------------------------------------------
# Function: Peeling algorithm for differentially private top-s selection
# --------------------------------------------
# Inputs:
#   m        - vector
#   sig      - sensitivity parameter
#   s        - number of top entries to select
#   epsilon  - privacy budget
#   delta    - privacy parameter (for approximate DP)
#   k        - dimension or a scaling parameter used in noise calibration
# Output:
#   p.set    - indices of the selected top-s entries under DP
peeling <- function(m , sig, s, epsilon,delta,k){
  nois <- rlaplace(n = m, location = 0, scale = sig*sqrt(3*k*s*log(1/delta))/epsilon)
  p.set<-c(which.max(m+nois))
  m[p.set] <- NA
  for (i in 1:(s-1)) {
    nois <- rlaplace(n = m, location = 0, scale = sig*sqrt(3*k*s*log(1/delta))/epsilon)
    p.set<-c(p.set,which.max((m+nois)))
    m[p.set]<-NA
  }
  return(p.set)
}
