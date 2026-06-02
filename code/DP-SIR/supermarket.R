# --- Load required libraries ---
library(LassoSIR)
library(glmnet)
library(geigen)
library(VGAM)
library(Matrix)
library(ROCR)
library(plotROC)
library(ggplot2)
library(pROC)
library(mgcv)
library(cowplot)
# Load all functions and variables defined in functions.R into the current R session
source("functions.R")
# --- Load dataset ---
market <- read.csv("market.csv", header=FALSE)
n <- dim(market)[1]   # number of observations
p <- dim(market)[2]   # number of variables
Y <- market[,1]       # response variable
X <- as.matrix(market[,-1])  # predictor matrix

n <- dim(X)[1]
p <- dim(X)[2]
Sigma.X <- cov(X)
Matrix.M <- matrix(data = 0, nrow = p, ncol = p)
ind <- as.numeric(Y<0)
quantile(Y,probs = c(0.33,0.66))
ind <- as.numeric(Y< -0.4424800)+as.numeric(Y<0.4148836)
# Feature screening
for (i in c(0,1,2)) {
  slicemean <- colMeans(X[which(ind==i),])
  Matrix.M  <- Matrix.M + (slicemean) %*% t(slicemean) * (length(which(ind==i))/n)
}
set.seed(951004)
dt_index <- order(diag(Matrix.M),decreasing = T)[1:20]
lasso_sir <- LassoSIR(X = X[,dt_index], Y = Y,H = 3, choosing.d = "automatic",categorical = FALSE,no.dim = 0,screening = TRUE)
lasso.sir.X <- as.vector(t(lasso_sir$beta)%*%t(X[,dt_index]))
lasso.sir.model <- mgcv::gam(Y~s(lasso.sir.X)) #(Table 3 in the main document)
summary(lasso.sir.model)

#####initial#####
scale_column <- function(x) {
  (2 * (x - min(x)) / (max(x) - min(x))) - 1
}

# Apply the transformation to each column
X_scaled <- X
X_scaled[which(X > 2)] <- 2
X_scaled[which(X < -2)] <- -2
K <- 1
H <- 7
s <- 7
# Quantize Y using DP histogram
slicey <- dp_hist(Y = Y, H = 5, epsilon = 0.1, m = 20)
set.seed(19951004)
epsilon <- 1
delta <- 1/n^{1}
tildey <- rep(1,n)
for (i in 1:(H-1)) {
  tildey[which(Y<=slicey[i])] <- tildey[which(Y<=slicey[i])] + 1
}
# Construct SIR matrix with scaled X
Sigma.X <- cov(X_scaled)
Martix.M <- matrix(data = 0, nrow = p, ncol = p)
overmean <- colMeans(X_scaled)
for (i in unique(tildey)[order(unique(tildey),decreasing = F)]) {
  if (length(which(tildey==i))==1) {slicemean <- X_scaled[which(tildey==i),]} else {
    slicemean <- colMeans(X_scaled[which(tildey==i),])}
  Martix.M  <- Martix.M + (slicemean-overmean) %*% t(slicemean-overmean) * (length(which(tildey==i))/n)
}
set.seed(1004)
# Peeling selection
pindex <- peeling(m = diag(Martix.M), sig = 6.01*2^2/n, s = s,epsilon = 15, delta = n^(-1.1),k=1)
pindex
# Add noise to compute private generalized eigen decomposition
set.seed(10093)
epsilon <- 2
delta <- n^(-1.1)
sigma.e1 <- 2*2*2^2*s/n*sqrt(2*log(2.5/delta))/epsilon
E1 <- rnorm(n = (s+1)*s/2, mean = 0, sd = sigma.e1)
E1_noise <- matrix(data = 0, nrow = s, ncol = s)
E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
E1_noise <- E1_noise+t(E1_noise)
diag(E1_noise)<- diag(E1_noise)/2
sigma.e2 <- 2*2*2^2*s/n*sqrt(2*log(2.5/delta))/epsilon
E2 <- rnorm(n = (s+1)*s/2, mean = 0, sd = sigma.e2)
E2_noise <- matrix(data = 0, nrow = s, ncol = s)
E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E1
E2_noise <- E2_noise+t(E2_noise)
diag(E2_noise)<- diag(E2_noise)/2
Sigma.X.noise.PD <- as.matrix(nearPD(Sigma.X[pindex,pindex] + E1_noise)$mat)
# Reconstruct full projection vector
ini.sir.estimator <- geigen(A = Martix.M[pindex,pindex] + E2_noise, B = Sigma.X.noise.PD, symmetric = T)
ini.sir.number <- order(ini.sir.estimator$values,decreasing = T)[1:K]
ini.sir.vector <- ini.sir.estimator$vectors[,ini.sir.number]
ini.sir.vector.complete <-matrix(data = 0, nrow = p, ncol = K)
ini.sir.vector.complete[pindex,] <- ini.sir.vector
ini.sir.X <- as.vector(t(ini.sir.vector.complete)%*%t(X))
# Fit GAM on projected covariates
ini.sir.model <- mgcv::gam(Y~s(ini.sir.X)) #(Table 3 in the main document)
summary(ini.sir.model)#0.72

# --- Estimate intrinsic dimension via DP-BIC ---
order.lambda <- ini.sir.estimator$values[order(ini.sir.estimator$values,decreasing = T)]
K = dp_bic(order.lambda = order.lambda, Cn = 2, HH = length(unique(tildey)))
K
#1

##dp.rf
set.seed(19951004)
TT<-2# Number of refinement iterations
dp.rf <- ini.sir.vector.complete# Initialize with previous DP-SIR estimate
dim(dp.rf) <- c(1,p)# Reshape for matrix multiplication
for (T.iter in 1:TT) {
  # Compute noise scale parameters for gradient and covariance perturbation
  sigma.e1 <- 2*2*4*s/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
  sigma.e2 <- 2*6*4*s/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
  # Generate symmetric noise matrices E1 and E2 for the signal and covariance matrices
  E1 <- rnorm(n = (p+1)*p/2, mean = 0, sd = sigma.e1)
  E1_noise <- matrix(data = 0, nrow = p, ncol = p)
  E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
  E1_noise <- E1_noise+t(E1_noise)
  diag(E1_noise)<- diag(E1_noise)
  E2 <- rnorm(n = (p+1)*p/2, mean = 0, sd = sigma.e2)
  E2_noise <- matrix(data = 0, nrow = p, ncol = p)
  E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E2
  E2_noise <- E2_noise+t(E2_noise)
  diag(E2_noise)<- diag(E2_noise)
  # Perturbed covariance matrix
  Sigma.X.noise.PD <- Sigma.X + E2_noise
  # Compute generalized Rayleigh quotient
  rhot <- as.numeric((dp.rf)%*%(Martix.M + E1_noise)%*%t(dp.rf)/((dp.rf)%*%as.matrix( Sigma.X.noise.PD)%*%t(dp.rf)))
  # Construct update matrix Ct for projection refinement
  Ct <- diag(p)+((Martix.M + E1_noise)-rhot* Sigma.X.noise.PD)/rhot/100
  # Update direction and project onto ℓ₀-constrained set
  dp.rf <- dp.rf%*%Ct/norm(dp.rf%*%Ct,type = "F")
  dp.rf[-order(abs(dp.rf),decreasing = T)[1:s]]<-0# Retain top-s entries
}
# Project data onto refined direction and fit GAM
dt.sir.X <- X%*%t(dp.rf)
dt.sir.model <- mgcv::gam(Y~s(dt.sir.X)) #(Table 3 in the main document)
summary(dt.sir.model)#0.715
###sir_prop####
set.seed(951004)
dp.sir <- ini.sir.vector.complete# Start from initial DP estimate
dim(dp.sir) <- c(p,1)
TT <- 1# Number of iterations
lambda <- 1# Regularization parameter
stepp <- 100# Learning rate denominator
for (T.iter in 1:TT) {
  # Compute noise scale for the gradient perturbation
  sigma.e1 <- (6*2+lambda*2*0.5)*3/n/stepp
  # Compute regularized gradient
  gd <- -Martix.M%*%dp.sir+lambda*Sigma.X%*%dp.sir%*%(t(dp.sir)%*%Sigma.X%*%dp.sir-diag(1))
  gd <- gd
  # Generate gradient noise
  E1 <- rnorm(n = s*1, mean = 0, sd = sigma.e1*sqrt(s)*sqrt(2*log(1.25*TT*6/5/delta))/epsilon*TT)
  # Gradient step
  dp.sir <- dp.sir - gd/stepp
  # Apply private feature selection using peeling
  dp.sir.index <- peeling(m = sqrt(rowSums(dp.sir^2)) , sig = sigma.e1, s = s, epsilon = epsilon/TT/6,delta = delta/TT/6,k=K)
  dp.sir.index
  # Inject noise and zero out unselected coefficients
  dp.sir[dp.sir.index,] <- dp.sir[dp.sir.index,] + E1
  dp.sir[-dp.sir.index,] <- 0
  
  # Fit model and evaluate deviance explained
  dp.sir.X <- X%*%(dp.sir)
  dp.sir.model <- mgcv::gam(Y~s(dp.sir.X))
  a <- summary(dp.sir.model)
  print(a$dev.expl)
}
# Final model after optimization
dp.sir.X <- X%*%(dp.sir)
dp.sir.model <- mgcv::gam(Y~s(dp.sir.X)) #(Table 3 in the main document)

which(dp.sir!=0)
summary(dp.sir.model)

# --- Plot fitted GAM with confidence intervals --- (Figure 1 in the main document)
p <- predict(dp.sir.model,se.fit = TRUE)
upr <- p$fit + p$se.fit*1.96
lor <- p$fit - p$se.fit*1.96
dat <- data.frame(x = dp.sir.X, y = Y,ex = dp.sir.model$fitted.values,up = upr, lo = lor) 
ggplot(data = dat)+geom_point(aes(x = x, y = y))+geom_line(aes(x = x, y = ex))+xlab("Covariates projected on the estimated differentially private direction")+ylab("Y")+
  geom_ribbon(aes(x = x, ymin=lor,ymax=upr),alpha=0.3)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))

# --- Plot marginal effect of X5 and X23 for comparison --- (Figure 2 in the main document)
x5 <- mgcv::gam(Y~s(X[,5]))
x23 <- mgcv::gam(Y~s(X[,23]))

# Plot for X5
examx <- x5
p <- predict(examx,se.fit = TRUE)
upr <- p$fit + p$se.fit*1.96
lor <- p$fit - p$se.fit*1.96
dat5 <- data.frame(x = X[,5], y = Y,ex = p$fit,up = upr, lo = lor) 
p5<-ggplot(data = dat5)+geom_point(aes(x = x, y = y))+geom_line(aes(x = x, y = ex))+xlab("X5")+ylab("Y")+
  geom_ribbon(aes(x = x, ymin=lo,ymax=up),alpha=0.3)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Plot for X23
examx <- x23
p <- predict(examx,se.fit = TRUE)
upr <- p$fit + p$se.fit*1.96
lor <- p$fit - p$se.fit*1.96
dat23 <- data.frame(x = X[,23], y = Y,ex = p$fit,up = upr, lo = lor) 
p23<-ggplot(data = dat23)+geom_point(aes(x = x, y = y))+geom_line(aes(x = x, y = ex))+xlab("X23")+ylab("Y")+
  geom_ribbon(aes(x = X[,23], ymin=lo,ymax=up),alpha=0.3)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Combine plots for visual comparison
plot_grid(p5, p23, nrow = 1)
