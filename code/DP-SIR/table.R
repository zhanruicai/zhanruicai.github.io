# --------------------------------------------
# Load required package
# --------------------------------------------
library(xtable)
# --------------------------------------------
# Load Low-Dimensional Simulation Results
# --------------------------------------------

# Linear Model
ld_linear_n2w_p15 <- read.csv("result/ld_linear_n2w_p15.csv")[,-1]
ld_linear_n2w_p30 <- read.csv("result/ld_linear_n2w_p30.csv")[,-1]
ld_linear_n4w_p15 <- read.csv("result/ld_linear_n4w_p15.csv")[,-1]
ld_linear_n4w_p30 <- read.csv("result/ld_linear_n4w_p30.csv")[,-1]

# Exponential Model
ld_exp_n2w_p15 <- read.csv("result/ld_exp_n2w_p15.csv")[,-1]
ld_exp_n2w_p30 <- read.csv("result/ld_exp_n2w_p30.csv")[,-1]
ld_exp_n4w_p15 <- read.csv("result/ld_exp_n4w_p15.csv")[,-1]
ld_exp_n4w_p30 <- read.csv("result/ld_exp_n4w_p30.csv")[,-1]

# M3 Model
ld_t1_n3w_p10 <- read.csv("result/ld_t1_n3w_p10.csv")[,-1]
ld_t1_n3w_p15 <- read.csv("result/ld_t1_n3w_p15.csv")[,-1]
ld_t1_n5w_p10 <- read.csv("result/ld_t1_n5w_p10.csv")[,-1]
ld_t1_n5w_p15 <- read.csv("result/ld_t1_n5w_p15.csv")[,-1]

# M4 Model
ld_t2_n3w_p10 <- read.csv("result/ld_t2_n3w_p10.csv")[,-1]
ld_t2_n3w_p15 <- read.csv("result/ld_t2_n3w_p15.csv")[,-1]
ld_t2_n5w_p10 <- read.csv("result/ld_t2_n5w_p10.csv")[,-1]
ld_t2_n5w_p15 <- read.csv("result/ld_t2_n5w_p15.csv")[,-1]

# Combine results into matrices by model
M1 <- rbind(rowMeans(ld_linear_n2w_p15),rowMeans(ld_linear_n2w_p30),rowMeans(ld_linear_n4w_p15),rowMeans(ld_linear_n4w_p30))
M2 <- rbind(rowMeans(ld_exp_n2w_p15),rowMeans(ld_exp_n2w_p30),rowMeans(ld_exp_n4w_p15),rowMeans(ld_exp_n4w_p30))
M3 <- rbind(rowMeans(ld_t1_n3w_p10),rowMeans(ld_t1_n3w_p15),rowMeans(ld_t1_n5w_p10),rowMeans(ld_t1_n5w_p15))
M4 <- rbind(rowMeans(ld_t2_n3w_p10),rowMeans(ld_t2_n3w_p15),rowMeans(ld_t2_n5w_p10),rowMeans(ld_t2_n5w_p15))

data <- rbind(M1,M2,M3,M4)

# Combine and display the full low-dimensional table
xtable(data,digits = 3)

# --------------------------------------------
# Load High-Dimensional Simulation Results
# --------------------------------------------

# Linear Model
hd_linear_n1k_p1k <- read.csv("result/hd_linear_n1k_p1k.csv")[,-1]
hd_linear_n1k_p2k <- read.csv("result/hd_linear_n1k_p2k.csv")[,-1]
hd_linear_n2k_p1k <- read.csv("result/hd_linear_n2k_p1k.csv")[,-1]
hd_linear_n2k_p2k <- read.csv("result/hd_linear_n2k_p2k.csv")[,-1]

# Exponential Model
hd_exp_n1k_p1k <- read.csv("result/hd_exp_n1k_p1k.csv")[,-1]
hd_exp_n1k_p2k <- read.csv("result/hd_exp_n1k_p2k.csv")[,-1]
hd_exp_n2k_p1k <- read.csv("result/hd_exp_n2k_p1k.csv")[,-1]
hd_exp_n2k_p2k <- read.csv("result/hd_exp_n2k_p2k.csv")[,-1]

# M3 Model
hd_t1_n2k_p2k <- read.csv("result/hd_t1_n2k_p2k.csv")[,-1]
hd_t1_n2k_p4k <- read.csv("result/hd_t1_n2k_p4k.csv")[,-1]
hd_t1_n4k_p2k <- read.csv("result/hd_t1_n4k_p2k.csv")[,-1]
hd_t1_n4k_p4k <- read.csv("result/hd_t1_n4k_p4k.csv")[,-1]

# M4 Model
hd_t2_n2k_p2k <- read.csv("result/hd_t2_n2k_p2k.csv")[,-1]
hd_t2_n2k_p4k <- read.csv("result/hd_t2_n2k_p4k.csv")[,-1]
hd_t2_n4k_p2k <- read.csv("result/hd_t2_n4k_p2k.csv")[,-1]
hd_t2_n4k_p4k <- read.csv("result/hd_t2_n4k_p4k.csv")[,-1]

# Combine results into matrices by model
M1 <- rbind(rowMeans(hd_linear_n1k_p1k),rowMeans(hd_linear_n1k_p2k),rowMeans(hd_linear_n2k_p1k),rowMeans(hd_linear_n2k_p2k))
M2 <- rbind(rowMeans(hd_exp_n1k_p1k),rowMeans(hd_exp_n1k_p2k),rowMeans(hd_exp_n2k_p1k),rowMeans(hd_exp_n2k_p2k))
M3 <- rbind(rowMeans(hd_t1_n2k_p2k,na.rm = T),rowMeans(hd_t1_n2k_p4k,na.rm = T),rowMeans(hd_t1_n4k_p2k,na.rm = T),rowMeans(hd_t1_n4k_p4k,na.rm = T))
M4 <- rbind(rowMeans(hd_t2_n2k_p2k,na.rm = T),rowMeans(hd_t2_n2k_p4k,na.rm = T),rowMeans(hd_t2_n4k_p2k,na.rm = T),rowMeans(hd_t2_n4k_p4k,na.rm = T))

# Combine and display the full high-dimensional table
data <- rbind(M1,M2,M3,M4)
xtable(data,digits = 3)
