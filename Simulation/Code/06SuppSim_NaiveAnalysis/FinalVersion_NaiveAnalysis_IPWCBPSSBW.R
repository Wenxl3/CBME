###############################################
# Simulation to Compare the Impact of Covariate Measurement Error
# on ATT Estimation using IPW, CBPS, SBW (ebcw) and SBW methods
#
# Comparison indices include:
#   - ATT bias and MSE,
#   - Balance metrics: ASMD (for X1 and U1) and WSMD.
#
# Note:
#   - Within each replication, weights are normalized within each group,
#     i.e., the treated group weights sum to 1 and the control group weights sum to 1.
#   - For balance assessment, the true (error-free) covariates are used.
#   - SBW is implemented via WeightIt with method = "ebcw", and SBW via sbw.
###############################################

# ------------------------- Environment Setup -------------------------
rm(list = ls())
print(getwd())
options(warn = -1)

# ------------------------- Package Loading -------------------------
library(MASS)       # For mvrnorm()
library(CBPS)       # For CBPS method
library(WeightIt)   # For weightit() - stable balancing weights by setting method = "optweight" in the SBWl to weightit().
library(openxlsx)   # For writing Excel files
# (ggplot2, ggpubr, reshape2 are not必需 for table output)

# ------------------------- Function Definitions -------------------------

# IPW Estimation Function with normalized weights
estimate_IPW <- function(data) {
  # Fit logistic regression using observed covariates: obsX and U
  model <- glm(W ~ obsX + U, data = data, family = binomial(link = "logit"))
  pi_hat <- predict(model, type = "response")

  # Raw weights: treated weight = 1; control: pi_hat/(1-pi_hat)
  Weights <- rep(NA, nrow(data))
  Weights[data$W == 1] <- 1
  Weights[data$W == 0] <- pi_hat[data$W == 0] / (1 - pi_hat[data$W == 0])

  # Normalize weights within each group (so that sum of weights equals 1)
  treated_idx <- which(data$W == 1)
  control_idx <- which(data$W == 0)
  Weights[treated_idx] <- Weights[treated_idx] / sum(Weights[treated_idx])
  Weights[control_idx] <- Weights[control_idx] / sum(Weights[control_idx])

  # ATT: difference between treated mean outcome and weighted control mean outcome
  ATT <- mean(data$Y[treated_idx]) - sum(Weights[control_idx] * data$Y[control_idx])
  return(list(ATT = ATT, Weights = Weights))
}

# CBPS Estimation Function with normalized weights
estimate_CBPS <- function(data) {
  # Fit CBPS model with ATT=1 so that treated weight is set to 1 by default
  cbps_fit <- CBPS(W ~ obsX + U, data = data,
                   # Set to 1 to find the average treatment effect on the treated interpreting the second level of the treatment factor as the treatment. .
                   ATT = 1,
                   method = "over",
                   standardize = TRUE)
  Weights <- cbps_fit$weights

  # # Standardize weights within each group
  treated_idx <- which(data$W == 1)
  control_idx <- which(data$W == 0)
  # Weights[treated_idx] <- Weights[treated_idx] / sum(Weights[treated_idx])
  # Weights[control_idx] <- Weights[control_idx] / sum(Weights[control_idx])

  ATT <- mean(data$Y[treated_idx]) - sum(Weights[control_idx] * data$Y[control_idx])
  return(list(ATT = ATT, Weights = Weights))
}

# SBW Estimation Function using WeightIt (method = "optweight") with normalized weights
## estimating optimization-based weights (also known as stable balancing weights)
estimate_SBW <- function(data) {
  SBW_fit <- weightit(W ~ obsX + U,
                      data = data,
                      method = "optweight", estimand = "ATT",
                      tols = 0)
  Weights <- SBW_fit$weights

  # Normalize weights within each group
  treated_idx <- which(data$W == 1)
  control_idx <- which(data$W == 0)
  Weights[treated_idx] <- Weights[treated_idx] / sum(Weights[treated_idx])
  Weights[control_idx] <- Weights[control_idx] / sum(Weights[control_idx])

  # ATT estimate
  ATT <- mean(data$Y[treated_idx]) - sum(Weights[control_idx] * data$Y[control_idx])
  return(list(ATT = ATT, Weights = Weights))
}

# ------------------------- Balance Metrics Functions -------------------------
# ASMD SBWculation Function
# SBWculates Absolute Standardized Mean Difference (ASMD) for each covariate
# using weighted group sums; uses true covariates (without measurement error).
# Global variables required:
#   - Z_NoME: matrix of true covariate values (columns: X and U)
#   - Index_trt: indices for treated subjects
#   - Index_col: indices for control subjects
#   - Z1_NoME: matrix of treated subjects (used to compute SD for normalization)
ASMDSBW <- function(Weights) {
  WeightedCovariate <- Z_NoME * Weights
  Num_trt <- colSums(WeightedCovariate[Index_trt,])
  Num_col <- colSums(WeightedCovariate[Index_col,])
  Denom <- apply(Z1_NoME, 2, sd)
  ASMD_values <- abs(Num_trt - Num_col) / Denom
  return(data.frame(X1ASMD = ASMD_values[1], U1ASMD = ASMD_values[2]))
}

# WSMD SBWculation Function
# Computes the Weighted Sample Mahalanobis Distance (WSMD) based on the difference
# of the weighted means of treated and control groups.
# Global variables:
#   - Z0_NoME: matrix of control group true covariates
#   - Z1_NoME: matrix of treated group true covariates
#   - Index_col: indices for control subjects
WSMDSBW <- function(Weights) {
  Y0_weighted <- Z0_NoME * Weights[Index_col]
  Mu1 <- colMeans(Z1_NoME)
  Mu0 <- colSums(Y0_weighted)
  Diff <- matrix(Mu1 - Mu0, ncol = 1)
  WSMD_val <- sqrt(t(Diff) %*% solve(var(Z1_NoME)) %*% Diff)
  return(as.numeric(WSMD_val))
}

# ------------------------- Basic Information -------------------------
p_Z <- 2                # Total number of covariates
p_X <- 1                # Number of error-prone covariates
p_U <- p_Z - p_X        # Number of error-free covariates
Mean_X <- 5             # True mean for error-prone covariate
Mean_U <- 10            # True mean for error-free covariate
Mu_Z <- c(Mean_X, Mean_U)
sd_X <- 1
sd_U1 <- 1

# Treatment assignment model parameters (for logistic model)
True.theta <- c(0.5, -3, 1.5)

# Outcome generation parameters:
# Model: Y = 210 + 27.4*X + 13.7*U + 10*W + error, where error ~ N(0, 4)
Meancoef <- c(210, 27.4, 13.7)
true_ATT <- 10
# Error variance (SD = 2)
Sigma_Y <- 4

# ------------------------- Simulation Information -------------------------
set.seed(102)
# For demonstration purposes，这里使用较小的 N 和 N_sim；实际应用中可调为 50000 和 1000
N <- 1000          # Sample size per simulation (adjust to 50000 as needed)
N_sim <- 1000       # Number of simulation replications (adjust to 1000 as needed)
rhoZ <- 0.3       # Correlation between covariates
vec_Vareps <- seq(0, 0.5, by = 0.1)  # Measurement error variances for error-prone covariate X

# ------------------------- Result Storage -------------------------
# 创建一个数据框存储每个测量误差水平和方法下的 ATT bias, MSE, ASMD及WSMD对比结果
results <- data.frame(
  Var_eps = rep(vec_Vareps, each = 3),
  Method = rep(c("IPW", "CBPS", "SBW"), times = length(vec_Vareps)),
  bias = NA,
  mse = NA,
  ASMD_X1 = NA,
  ASMD_U1 = NA,
  WSMD = NA
)

# 如有需要，也可保存每个测量误差水平下的详细结果
result_list <- list()

# ------------------------- Main Simulation Loop -------------------------
for (ve in vec_Vareps) {
  # 初始化用于存储当前测量误差水平下各方法结果的向量/矩阵
  ATT_IPW_vec <- numeric(N_sim)
  ATT_CBPS_vec <- numeric(N_sim)
  ATT_SBW_vec <- numeric(N_sim)

  ASMD_IPW_mat <- matrix(0, nrow = N_sim, ncol = 2)  # Columns: X1ASMD, U1ASMD
  WSMD_IPW_vec <- numeric(N_sim)

  ASMD_CBPS_mat <- matrix(0, nrow = N_sim, ncol = 2)
  WSMD_CBPS_vec <- numeric(N_sim)

  ASMD_SBW_mat <- matrix(0, nrow = N_sim, ncol = 2)
  WSMD_SBW_vec <- numeric(N_sim)

  for (i in 1:N_sim) {
    # -------------------- Data Generation --------------------
    # Generate true covariate matrix Z ~ N(Mu_Z, Sigma_Z)
    Sigma_Z <- matrix(c(sd_X^2, rhoZ * sd_X * sd_U1,
                        rhoZ * sd_X * sd_U1, sd_U1^2),
                      ncol = p_Z, byrow = TRUE)
    Z <- mvrnorm(n = N, mu = Mu_Z, Sigma = Sigma_Z)
    colnames(Z) <- c("X", "U")

    # Generate measurement error for X and produce observed covariate obsX
    eps_X <- rnorm(n = N, mean = 0, sd = sqrt(ve))
    obsX <- Z[, "X"] + eps_X  # observed (error-prone) X

    # Treatment assignment: generate propensity scores using true covariates and assign treatment
    eta <- as.numeric(Z %*% True.theta[-1] + True.theta[1])
    pi_true <- exp(eta) / (1 + exp(eta))
    W <- rbinom(n = N, size = 1, prob = pi_true)

    # Outcome generation: potential outcomes and observed outcome Y
    Y0 <- Meancoef[1] +
      Meancoef[2] * Z[, "X"] +
      Meancoef[3] * Z[, "U"]
    Y1 <- Y0 + true_ATT
    error_term <- rnorm(n = N, mean = 0, sd = sqrt(Sigma_Y))
    Y <- ifelse(W == 1, Y1, Y0) + error_term

    # Assemble observed data
    data <- data.frame(obsX = obsX, U = Z[, "U"], W = W, Y = Y)

    # -------------------- ATT Estimation --------------------
    # Apply each method to obtain ATT estimates and corresponding weights
    res_ipw <- estimate_IPW(data)
    ATT_IPW_vec[i] <- res_ipw$ATT
    weights_ipw <- res_ipw$Weights

    res_cbps <- estimate_CBPS(data)
    ATT_CBPS_vec[i] <- res_cbps$ATT
    weights_cbps <- res_cbps$Weights

    res_SBW <- estimate_SBW(data)  # SBW implemented via WeightIt (method = "ebcw")
    ATT_SBW_vec[i] <- res_SBW$ATT
    weights_SBW <- res_SBW$Weights

    # -------------------- Balance Metrics SBWculation --------------------
    # For balance metrics，use the true covariate matrix (without measurement error)
    Z_NoME <- as.matrix(Z)
    Index_trt <- which(W == 1)
    Index_col <- which(W == 0)
    Z1_NoME <- Z_NoME[Index_trt, , drop = FALSE]
    Z0_NoME <- Z_NoME[Index_col, , drop = FALSE]

    # IPW balance metrics
    asmd_ipw <- ASMDSBW(weights_ipw)
    wsmd_ipw <- WSMDSBW(weights_ipw)
    ASMD_IPW_mat[i,] <- c(asmd_ipw$X1ASMD, asmd_ipw$U1ASMD)
    WSMD_IPW_vec[i] <- wsmd_ipw

    # CBPS balance metrics
    asmd_cbps <- ASMDSBW(weights_cbps)
    wsmd_cbps <- WSMDSBW(weights_cbps)
    ASMD_CBPS_mat[i,] <- c(asmd_cbps$X1ASMD, asmd_cbps$U1ASMD)
    WSMD_CBPS_vec[i] <- wsmd_cbps

    # SBW balance metrics
    asmd_SBW <- ASMDSBW(weights_SBW)
    wsmd_SBW <- WSMDSBW(weights_SBW)
    ASMD_SBW_mat[i,] <- c(asmd_SBW$X1ASMD, asmd_SBW$U1ASMD)
    WSMD_SBW_vec[i] <- wsmd_SBW

    print(paste0("ve = ", ve, ", ith-Sim = ", i))
  }  # End simulation replications loop

  # Compute summary statistics for current measurement error variance (ve)
  bias_ipw <- mean(ATT_IPW_vec - true_ATT)
  mse_ipw <- mean((ATT_IPW_vec - true_ATT)^2)
  asmd_ipw_mean <- colMeans(ASMD_IPW_mat)
  wsmd_ipw_mean <- mean(WSMD_IPW_vec)

  bias_cbps <- mean(ATT_CBPS_vec - true_ATT)
  mse_cbps <- mean((ATT_CBPS_vec - true_ATT)^2)
  asmd_cbps_mean <- colMeans(ASMD_CBPS_mat)
  wsmd_cbps_mean <- mean(WSMD_CBPS_vec)

  bias_SBW <- mean(ATT_SBW_vec - true_ATT)
  mse_SBW <- mean((ATT_SBW_vec - true_ATT)^2)
  asmd_SBW_mean <- colMeans(ASMD_SBW_mat)
  wsmd_SBW_mean <- mean(WSMD_SBW_vec)

  # 保存当前测量误差水平下的详细结果到列表中（如有需要）
  result_list[[as.character(ve)]] <- list(
    ATT_IPW = ATT_IPW_vec, bias_IPW = bias_ipw, mse_IPW = mse_ipw,
    ASMD_IPW = asmd_ipw_mean, WSMD_IPW = wsmd_ipw_mean,

    ATT_CBPS = ATT_CBPS_vec, bias_CBPS = bias_cbps, mse_CBPS = mse_cbps,
    ASMD_CBPS = asmd_cbps_mean, WSMD_CBPS = wsmd_cbps_mean,

    ATT_SBW = ATT_SBW_vec, bias_SBW = bias_SBW, mse_SBW = mse_SBW,
    ASMD_SBW = asmd_SBW_mean, WSMD_SBW = wsmd_SBW_mean
  )

  # 将当前测量误差水平下的统计摘要存入 results 数据框中
  results$bias[results$Var_eps == ve & results$Method == "IPW"] <- bias_ipw
  results$mse[results$Var_eps == ve & results$Method == "IPW"] <- mse_ipw
  results$ASMD_X1[results$Var_eps == ve & results$Method == "IPW"] <- asmd_ipw_mean[1]
  results$ASMD_U1[results$Var_eps == ve & results$Method == "IPW"] <- asmd_ipw_mean[2]
  results$WSMD[results$Var_eps == ve & results$Method == "IPW"] <- wsmd_ipw_mean

  results$bias[results$Var_eps == ve & results$Method == "CBPS"] <- bias_cbps
  results$mse[results$Var_eps == ve & results$Method == "CBPS"] <- mse_cbps
  results$ASMD_X1[results$Var_eps == ve & results$Method == "CBPS"] <- asmd_cbps_mean[1]
  results$ASMD_U1[results$Var_eps == ve & results$Method == "CBPS"] <- asmd_cbps_mean[2]
  results$WSMD[results$Var_eps == ve & results$Method == "CBPS"] <- wsmd_cbps_mean

  results$bias[results$Var_eps == ve & results$Method == "SBW"] <- bias_SBW
  results$mse[results$Var_eps == ve & results$Method == "SBW"] <- mse_SBW
  results$ASMD_X1[results$Var_eps == ve & results$Method == "SBW"] <- asmd_SBW_mean[1]
  results$ASMD_U1[results$Var_eps == ve & results$Method == "SBW"] <- asmd_SBW_mean[2]
  results$WSMD[results$Var_eps == ve & results$Method == "SBW"] <- wsmd_SBW_mean

  cat("Completed Var_eps =", ve, "\n")
}  # End outer loop over vec_Vareps

# ------------------------- Final Output -------------------------
print("Summary Comparison Table (by Measurement Error Variance and Method):")
print(results)

# 将 results 中所有数值型变量转换为格式化字符串（保留三位小数）
results_formatted <- results
numeric_cols <- sapply(results_formatted, is.numeric)
results_formatted[numeric_cols] <- lapply(results_formatted[numeric_cols],
                                          function(x) sprintf("%.3f", x))


# Export results as an Excel file
filename <- paste0("Supp_NaiveAnalysis_N", N / 1000, "k_Sim", N_sim / 1000, "k_IPW_CBPS_SBW.xlsx")
write.xlsx(results_formatted, file = filename, rowNames = FALSE)
