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

##### --- Peeling function for differentially private variable selection ---
peeling <- function(m , sig, s){
  nois <- rlaplace(n = m, location = 0, scale = sig)
  p.set<-c(which.max(m+nois))
  m[p.set] <- NA
  for (i in 1:(s-1)) {
    nois <- rlaplace(n = m, location = 0, scale = sig)
    p.set<-c(p.set,which.max((m+nois)))
    m[p.set]<-NA
  }
  return(p.set)
}
### --- Load Arcene dataset (UCI) ---
#https://archive.ics.uci.edu/ml/machine-learning-databases/arcene/
arcene_train <- read.table("arcene/ARCENE/arcene_train.data", quote="\"", comment.char="")
arcene_test <- read.table("arcene/ARCENE/arcene_test.data", quote="\"", comment.char="")
arcene_valid <- read.table("arcene/ARCENE/arcene_valid.data", quote="\"", comment.char="")
arcene_trainy <- read.table("arcene/ARCENE/arcene_train.labels", quote="\"", comment.char="")
arcene_valid <- read.table("arcene/ARCENE/arcene_valid.data", quote="\"", comment.char="")
arcene_validy <- read.table("arcene/arcene_valid.labels", quote="\"", comment.char="")

### --- Preprocessing ---
set.seed(951004)
X <- data.matrix(arcene_train, rownames.force = NA)
pX <- data.matrix(arcene_test, rownames.force = NA)
# Fit logistic Lasso on training set
lasso_glm.cv <- cv.glmnet(x = X, y = (arcene_trainy$V1+1)/2, family = "binomial")
lasso_glm <- glmnet(x = X, y = (arcene_trainy$V1+1)/2, family = "binomial",lambda = lasso_glm.cv$lambda.min)
line_pred <- lasso_glm$a0+as.numeric(lasso_glm$beta)%*%t(pX)
predictions_lasso_sir <- 1/(1+exp(-line_pred))
# Generate probabilistic predictions and sample pseudo-labels for test data
py <- vapply(predictions_lasso_sir, function(x) rbinom(1, 1, x), as.integer(1L))
# Combine training and pseudo-labeled test data
nX<-rbind(X,pX)
ny<-c((arcene_trainy$V1+1)/2,py)
##n=100,p=10000
# Standardize full design matrix
X.s <- scale(nX)
# Standardize validation data using training-test combined statistics
VX <- data.matrix(arcene_valid, rownames.force = NA)
VX.s <- VX%*%diag(1/sqrt(apply(nX, 2, sd)))
VY <- data.matrix(arcene_validy, rownames.force = NA)
### --- Benchmark GLM using logistic Lasso ---
set.seed(951004)
cv <- cv.glmnet(x = X.s, y = ny,family = "binomial" ,standardize = FALSE, type.measure = "auc")
glm_lasso <- glmnet(x = X.s,y = ny,family = "binomial",standardize=FALSE, lambda = cv$lambda.min)
### --- Lasso-SIR ---
set.seed(951004)
lasso_sir <- LassoSIR(X = X.s, Y = ny,H=0,categorical = TRUE,no.dim = 1)
lasso_sir_glm <- glm(ny~X.s%*%lasso_sir$beta,family = "binomial")
### --- SIR with differentially private initialization ---
# Truncate values for stability
X.s[which(X.s > 1)]<-1 
X.s[which(X.s < -1)]<--1 

n <- dim(X.s)[1]
p <- dim(X.s)[2]
# Compute between-slice covariance matrix
Sigma.X <- cov(X.s)
Matrix.M <- matrix(data = 0, nrow = p, ncol = p)
for (i in c(0,1)) {
  slicemean <- colMeans(X.s[which(ny==i),])
  Matrix.M  <- Matrix.M + (slicemean) %*% t(slicemean) * (length(which(ny==i))/n)
}
##### --- Initial private projection estimation using DP generalized eigenvalue problem ---
K <- 1   # Dimension of projection
H <- 2   # Number of slices (binary in this case)
s <- 5   # Sparsity level
set.seed(9510)
epsilon <- 2
delta <- 1/n^{1}
# Select 2s features using private peeling
pindex <- peeling(m = diag(Matrix.M), sig = 6.01/n*sqrt(2*s*log(1/delta))/epsilon, s = 2*s)
# Add Gaussian noise to Matrix.M and Sigma.X for privacy
sigma.e1 <- 2*2*s/n*sqrt(2*log(2.5/delta))/epsilon
E1 <- rnorm(n = (2*s+1)*2*s/2, mean = 0, sd = sigma.e1)
E1_noise <- matrix(data = 0, nrow = 2*s, ncol = 2*s)
E1_noise[lower.tri(E1_noise, diag=TRUE)] <- E1
E1_noise <- E1_noise+t(E1_noise)
diag(E1_noise)<- diag(E1_noise)/2
sigma.e2 <- 2*2*s/n*sqrt(2*log(2.5/delta))/epsilon
E2 <- rnorm(n = (2*s+1)*2*s/2, mean = 0, sd = sigma.e2)
E2_noise <- matrix(data = 0, nrow = 2*s, ncol = 2*s)
E2_noise[lower.tri(E2_noise, diag=TRUE)] <- E2
E2_noise <- E2_noise+t(E2_noise)
diag(E2_noise)<- diag(E2_noise)/2
# Ensure positive definiteness and solve generalized eigenvalue problem
Sigma.X.noise.PD <- as.matrix(nearPD(Sigma.X[pindex,pindex] + E1_noise)$mat)
# Extract top eigenvector and expand to full dimension
ini.sir.estimator <- geigen(A = Matrix.M[pindex,pindex] + E2_noise, B = Sigma.X.noise.PD, symmetric = T)
ini.sir.number <- order(ini.sir.estimator$values,decreasing = T)[1:K]
ini.sir.vector <- ini.sir.estimator$vectors[,ini.sir.number]
ini.sir.vector.complete <-matrix(data = 0, nrow = p, ncol = K)
ini.sir.vector.complete[pindex,] <- ini.sir.vector
# Fit logistic regression on projection
ini_sir_glm <- glm(ny~X.s%*%(ini.sir.vector.complete),family = "binomial")
# Predict on validation set using estimated direction
line_pred <- ini_sir_glm$coefficients[1]+ini_sir_glm$coefficients[-1]*as.vector(t(ini.sir.vector.complete)%*%t(VX.s))
predictions_ini_sir <- 1/(1+exp(-line_pred))
pred <- prediction(predictions_ini_sir,(VY+1)/2)
# Evaluate AUC using ROCR
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR2 <- auc_ROCR@y.values[[1]]
auc_ROCR2

##dp.rf
set.seed(951004)
TT <- 3  # Number of refinement iterations
dp.rf <- ini.sir.vector.complete  # Initialize with previous DP-SIR estimate
dim(dp.rf) <- c(1, p)  # Reshape as a row vector

for (T.iter in 1:TT) {
  # Compute privacy-preserving noise scales
  sigma.e1 <- 2*2*s/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
  sigma.e2 <- 2*6.01*s/n*sqrt(2*log(2.5*TT/delta))/epsilon*TT
  # Add symmetric Gaussian noise to signal and covariance matrices
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
  # Perturbed covariance matrix
  Sigma.X.noise.PD <- Sigma.X + E1_noise
  # Compute generalized Rayleigh quotient (signal-to-noise)
  rhot <- as.numeric((dp.rf)%*%(Matrix.M + E2_noise)%*%t(dp.rf)/((dp.rf)%*%as.matrix( Sigma.X.noise.PD)%*%t(dp.rf)))
  # Construct update matrix Ct using scaled eigenstructure
  Ct <- diag(p)+((Matrix.M + E2_noise)-rhot* Sigma.X.noise.PD)/rhot/2 ##step size:1/4
  # Gradient step and sparsification
  dp.rf <- dp.rf%*%Ct/norm(dp.rf%*%Ct,type = "F")
  dp.rf[-order(abs(dp.rf),decreasing = T)[1:s]]<-0
}
# Fit logistic regression on refined projection
dt_sir_glm <- glm(ny~X.s%*%t(dp.rf),family = "binomial")
# Predict on validation set and compute AUC
line_pred <- dt_sir_glm$coefficients[1]+dt_sir_glm$coefficients[2]*as.vector((VX.s)%*%t(dp.rf))
pred_dt_sir_glm <- 1/(1+exp(-line_pred))
pred <- prediction(pred_dt_sir_glm,(VY+1)/2)
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR3<- auc_ROCR@y.values[[1]]
auc_ROCR3
###sir_prop####
set.seed(1004)
dp.sir <- t(ini.sir.vector.complete)/norm(ini.sir.vector.complete)
dim(dp.sir) <- c(p,1)
lambda = 1
stepp <- 10
TT <- 2
for (T.iter in 1:TT) {
  # Compute gradient and noise scale
  sigma.e1 <- (6.01+lambda*2*1)/n*sqrt(K)/stepp
  gd <- -Matrix.M%*%dp.sir+lambda*Sigma.X%*%dp.sir%*%(t(dp.sir)%*%Sigma.X%*%dp.sir-diag(K))
  gd <- gd
  # Add noise to selected components
  E1 <- rnorm(n = s*K, mean = 0, sd = 6/5*sigma.e1*sqrt(s)*sqrt(2*log(1.25*TT*6/5/delta))/epsilon*TT)
  dp.sir <- dp.sir - gd/stepp
  # Perform private selection via peeling
  dp.sir.index <- peeling(m = sqrt(rowSums(dp.sir^2)) , sig = sigma.e1*sqrt(3*K*s*log(1/(delta/TT/6)))/(epsilon/TT/6), s = s)
  dp.sir.index
  dp.sir[dp.sir.index,] <- dp.sir[dp.sir.index,] + E1
  dp.sir[-dp.sir.index,] <- 0
  # Normalize to maintain bounded ℓ₂-norm
  for (i in 1:K) {
    if (norm(dp.sir[,i],type = "2")>2) {dp.sir[,i]<-dp.sir[,i]/norm(dp.sir[,i],type = "2")*2} 
  }
}
# Fit and evaluate logistic regression with DP-SIR
prop_sir_glm <- glm(ny~(X.s%*%dp.sir),family = "binomial")$coefficients
line_pred <- prop_sir_glm[1]+ prop_sir_glm[2]*as.vector((VX.s)%*%(dp.sir))
predictions_prop_sir_glm <- 1/(1+exp(-line_pred))
pred <- prediction(predictions_prop_sir_glm,(VY+1)/2)
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR4 <- auc_ROCR@y.values[[1]]
auc_ROCR4
### --- ROC Curve Comparison Across Methods ---
# Lasso-GLM

line_pred <- glm_lasso$a0+(VX.s)%*%data.matrix(glm_lasso$beta)
predictions_lasso_glm <- 1/(1+exp(-line_pred))
pred <- prediction(predictions_lasso_glm,(VY+1)/2)
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR1 <- auc_ROCR@y.values[[1]]
auc_ROCR1
rocobj1 <- roc(as.vector(VY+1)/2, as.vector(predictions_lasso_glm))
# Initial SIR
line_pred <- ini_sir_glm$coefficients[1]+ini_sir_glm$coefficients[-1]*as.vector(t(ini.sir.vector.complete)%*%t(VX.s))
predictions_ini_sir <- 1/(1+exp(-line_pred))
pred <- prediction(predictions_ini_sir,(VY+1)/2)
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR2 <- auc_ROCR@y.values[[1]]
auc_ROCR2
rocobj2 <- roc(as.vector(VY+1)/2, as.vector(predictions_ini_sir))
# DP-RF
line_pred <- dt_sir_glm$coefficients[1]+dt_sir_glm$coefficients[2]*as.vector((VX.s)%*%t(dp.rf))
pred_dt_sir_glm <- 1/(1+exp(-line_pred))
pred <- prediction(pred_dt_sir_glm,(VY+1)/2)
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR3<- auc_ROCR@y.values[[1]]
auc_ROCR3
rocobj3 <- roc(as.vector(VY+1)/2, as.vector(pred_dt_sir_glm))
# Proposed DP-SIR 
line_pred <- prop_sir_glm[1]+ prop_sir_glm[2]*as.vector((VX.s)%*%(dp.sir))
predictions_prop_sir_glm <- 1/(1+exp(-line_pred))
pred <- prediction(predictions_prop_sir_glm,(VY+1)/2)
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR4 <- auc_ROCR@y.values[[1]]
auc_ROCR4
rocobj4 <- roc(as.vector(VY+1)/2, as.vector(predictions_prop_sir_glm))

# Lasso-SIR
lasso_sir <- LassoSIR(X = X.s, Y = ny,H=0,categorical = TRUE,no.dim = 1)
lasso_sir_glm <- glm(ny~X.s%*%lasso_sir$beta,family = "binomial")
line_pred <- lasso_sir_glm$coefficients[1]+lasso_sir_glm$coefficients[-1]*as.vector(t(lasso_sir$beta)%*%t(VX.s))
predictions_lasso_sir <- 1/(1+exp(-line_pred))
pred <- prediction(predictions_lasso_sir,(VY+1)/2)
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR5 <- auc_ROCR@y.values[[1]]
auc_ROCR5
rocobj5 <- roc(as.vector(VY+1)/2, as.vector(predictions_lasso_sir))

# --- Create and plot ROC curves using ggplot2 ---
data.ggplot <- data.frame(Specificity = c(rocobj1$specificities,rocobj2$specificities,rocobj3$specificities,rocobj4$specificities,rocobj5$specificities),Sensitivity = c(rocobj1$sensitivities,rocobj2$sensitivities,rocobj3$sensitivities,rocobj4$sensitivities,rocobj5$sensitivities),method = c(rep("Lasso-glm",length(rocobj1$sensitivities)),rep("Lasso-SIR",length(rocobj2$sensitivities)),rep("DP-TR",length(rocobj3$sensitivities)),rep("Proposed DP-SIR",length(rocobj4$sensitivities)),rep("Lasso-SIR",length(rocobj5$sensitivities))))
ggplot(data = data.ggplot,aes(x=Specificity,y=Sensitivity))+geom_line(aes(color = method,linetype = method))+
  xlim(1,0)+ylim(0,1)
# Plot multiple ROC curves with AUC annotations
pROC::ggroc(list("DP-SIni; AUC:0.638" = rocobj2,"DP-TRF; AUC:0.684" = rocobj3, "DP-SSIR; AUC:0.773" = rocobj4, "Lasso-SIR; AUC:0.838" = rocobj5),aes = "color") + labs(color = 'Method')+ scale_colour_manual(values = c("red", "blue", "black","green"))+
  theme(legend.position = c(0.8, 0.2))+ geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed")
