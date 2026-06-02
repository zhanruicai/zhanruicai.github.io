# Differentially Private Sliced Inverse Regression: Minimax Optimality and Algorithm

## 1. Data

### Abstract

We use two real datasets in the analysis:
1. **The supermarket dataset** – Each record corresponds to a daily observation collected from a major supermarket in northern China. The response of interest is the number of customers on a particular day, and each predictor represents the sales volume of a specific product on that day.
2. **The Arcene dataset** – This dataset is used to distinguish cancer from normal patterns in mass-spectrometric data. It is a two-class classification problem with continuous input variables.

### Availability

The data files are not included in the supplementary materials. However, the Arcene dataset is publicly available and can be accessed online:

Guyon, I., Gunn, S., Ben-Hur, A., & Dror, G. (2004). *Arcene Dataset*. UCI Machine Learning Repository. Retrieved from [https://archive.ics.uci.edu/dataset/167/arcene](https://archive.ics.uci.edu/dataset/167/arcene).

The supermarket dataset cannot be made publicly available due to data sharing restrictions.

---

## 2. Code

### Abstract

All code is written in R and organized into modular scripts. It includes core functions for:
- Differentially private histogram slicing
- Sparse variable selection via noisy peeling
- Generalized eigenvalue decomposition with added noise
- Performing (sparse) differentially private SIR

### Availability

The full R code for this paper is provided, including:
- One file containing necessary functions
- Two files for the simulation settings in the main manuscript
- Two files for real data analyses
- Two additional scripts for DP histogram slicing and DP cross-validation (referenced in the supplementary materials)

### Description

The following R scripts are included for simulations, empirical analyses, and result visualization:

- **functions.R**
  Contains utility functions for:
  - Computing projection matrices
  - Estimating quantiles based on histogram CDFs
  - Differentially private histogram estimation
  - Generating autoregressive covariance matrices
  - Peeling algorithm for differentially private top‑s selection

- **low_dim.R**
  Generates simulation results for Table 1 in the main manuscript (high-dimensional settings). Results are stored in the `result` folder.

- **high_dim.R**
  Generates simulation results for Table 2 in the main manuscript (low-dimensional settings). Results are stored in the `result` folder.

- **table.R**
  Summarizes outcomes from `low_dim.R` and `high_dim.R` and compiles the final result tables for inclusion as Table 1 and Table 2 in the manuscript.

- **supermarket.R**
  Applies DP-SIR to the supermarket dataset, evaluating predictive performance using generalized additive models (GAMs). Supports Figure 1, Figure 2, and Table 3 in the main document.

- **arcene.R**
  Applies DP-SIR to the Arcene dataset for binary classification, benchmarking performance via logistic regression and ROC-AUC analysis. Supports Figure C.1 in the supplementary material.

- **illustration_hist.R**
  Illustrates the effects of differentially private histogram slicing on eigenstructure and subspace error. Supports Figures C.2–C.11 in the supplementary material.

- **tuns.R**
  Compares the proposed cross-validation estimator with an oracle benchmark. Supports Figure C.12 in the supplementary material.

---

## 3. Instructions for Use

### Reproducibility Scope

All plots and tables in the paper can be reproduced using the scripts above. Simulation outputs are stored in the `result` folder.

### Requirements

- **R version:** ≥ 4.0
- **Required packages:**
  ```r
  install.packages(c("geigen", "VGAM", "Matrix", "MASS", "glmnet", "ROCR", "pROC",
                     "plotROC", "mgcv", "cowplot", "foreach", "doParallel",
                     "xtable", "patchwork"))
  ```

Some scripts assume a registered parallel backend:
```r
cl_size <- 10
cl <- makeCluster(cl_size)
registerDoParallel(cl)
```
