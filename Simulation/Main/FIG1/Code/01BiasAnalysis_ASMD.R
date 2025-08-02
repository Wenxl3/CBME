# -*- coding: utf-8 -*-
# ------------------------------------------
# Title       : 01BiasAnalysis_ASMD.R
# Objective   : visualize the impact of measurement error in covariates on EB method
# Reference(s): Entropr balancing
# ------------------------------------------

# Calculate/Approximate the 【true】 cov(Z | T = 1)
{
{
  # Working Directory
  rm(list = ls())
  getwd()
  # Suppressing warnings in the global settings
  options(warn = -1)

  # Package
  library(MASS)

  # Basic Information----
  # covariates
  ## Dimension of all covariates
  p_Z <- 2
  ## Dimension of error-prone covariates
  p_X <- 1
  ## Dimension of error-free covariates
  p_U <- p_Z - p_X
  ## expectation of error-prone data
  Mean_X <- 5
  ## expectation of error-free data
  Mean_U <- 10
  ## Expectation vector of covariates
  Mu_Z <- c(Mean_X, Mean_U)
  sd_X <- 1
  sd_U1 <- 1
  # True logistic model parameter
  True.theta <- c(0.5, -3, 1.5)

  # Simulation Information ----
  set.seed(102)

  N <- 1000
  N_sim <- 1000
  rhoZ <- 0.3

  vec_Vareps <- seq(0, 0.5, by = 0.05)

  # Save Rdata: the estimate of Cov(Z|T=1)
  list_Cov_trt <- list()
}

  for (ii in 1:length(vec_Vareps)) {
    # ii: indicator of the i(th) Var(eps)
    ## ii <- 5
    Var_eps <- vec_Vareps[ii]

    SimSum_Var <- matrix(data = 0, nrow = p_Z, ncol = p_Z)

    for (i in 1:N_sim) {
      # i <- 1

      # Covariance matrix of covariates
      Sigma_Z <- matrix(c(sd_X^2, rhoZ * sd_X * sd_U1,
                          rhoZ * sd_X * sd_U1, sd_U1^2),
                        ncol = p_Z,
                        byrow = T)
      Sigma_eps <- diag(x = c(rep(Var_eps, p_X), rep(0, p_U)))

      # Covariates----
      ## ·Z----
      Z <- mvrnorm(n = N, mu = Mu_Z, Sigma = Sigma_Z)
      eps_X <- rnorm(n = N, mean = 0, sd = sqrt(Var_eps))
      colnames(Z) <- c(paste0("X", 1:p_X), paste0("U", 1:p_U))
      ### error-prone covariate
      X <- Z[, paste0("X", 1:p_X)]
      obsX <- X + eps_X
      ## ·Propensity score & # Treatment assignment----
      eta <- Z %*% True.theta[-1] + True.theta[1]
      pi <- exp(eta) / (1 + exp(eta))
      W <- sapply(X = pi, FUN = function(x) { rbinom(1, size = 1, x) })
      Index_trt <- which(W == 1)

      SimSum_Var <- SimSum_Var + var(Z[Index_trt,])

      print(paste0("Var(eps) = ", Var_eps, ", No.sim = ", i))
    }

    list_Cov_trt[[ii]] <- SimSum_Var / N_sim
  }

  # save(list_Cov_trt, file = "list_Cov_trt.RData")

  True_Cov_trt <- Reduce("+", list_Cov_trt) / length(vec_Vareps)
}

{

{
  # rm(list = ls())
  allvariable <- ls()
  rm(list = allvariable[which(allvariable != 'True_Cov_trt')])
  rm(allvariable)
  getwd()
  # Suppressing warnings in the global settings
  options(warn = -1)

  # Package Loading----
  library(MASS)
  library(ebal)
  # library(WeightIt)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)

  # Objective, Gradient, Hessian Functions----
  ## ·EB — L&W
  objective.zqy <- function(theta) {
    f <- log((1 - W) %*% exp(Z %*% theta)) - Z1_bar %*% theta
    return(f)
  }

  gradient.zqy <- function(theta) {
    w <- drop((1 - W) * exp(Z %*% theta))
    f <- drop(w %*% (Z) / sum(w) - Z1_bar)
    return(f)
  }

  hessian.zqy <- function(theta) {
    f <- 0
    w <- drop((1 - W) * exp(Z %*% theta))
    a <- matrix(rep(0, p_Z^2), nrow = p_Z, ncol = p_Z)
    for (j in 1:N) {
      a <- a + (Z[j,]) %*% t(Z[j,]) * w[j]
    }
    w2 <- a * sum(w) - t(w %*% Z) %*% (w %*% Z)
    f <- f + drop(w2 / (sum(w))^2)
    return(f)
  }

  # ASMD Calculation Function----
  ASMDCal <- function(Weights) {
    # Weights <- Weights_BA

    # Weighted Covariate (covariate without measurement error)
    WeightedCovariate <- Z_NoME * Weights

    # ASMD for each Covariate
    ## Numerator
    ASMD.Numerator_trt <- colSums(WeightedCovariate[Index_trt,])
    ASMD.Numerator_col <- colSums(WeightedCovariate[Index_col,])
    ## Denominator
    ASMD.Denominator1 <- apply(Z1_NoME, MARGIN = 2, sd)

    ##ASMD
    ASMD <- abs(ASMD.Numerator_trt - ASMD.Numerator_col) / ASMD.Denominator1

    return(data.frame(X1ASMD = ASMD[1], U1ASMD = ASMD[2]))
  }

  # Approximate ASMD Calculation Function----
  approxASMDCal <- function(Weights) {
    # Weights <- Weights_BA

    # Approximate ASMD for each Covariate
    ## Numerator
    approxASMD.Numerator <- Sigma_eps %*% par_BA
    ## Denominator
    ### Use true the variance
    approxASMD.Denominator1 <- sqrt(diag(True_Cov_trt))

    ##ASMD
    approxASMD <- abs(approxASMD.Numerator) / approxASMD.Denominator1

    return(data.frame(X1ASMD = approxASMD[1], U1ASMD = approxASMD[2]))
  }

  # Weighted Sample Mahalanobis distance(WSMD) Calculation Function----
  ## Using the function in "WSMDcal_UsingRInternalFunctions.R"
  WSMDCal <- function(Weights) {
    # Weights <- Weights_BA

    # Weighted covariate
    Y0.X_weighted <- Z0_NoME * Weights[Index_col]
    weighted_Mu1 <- colMeans(Z1_NoME)
    weighted_Mu0 <- colSums(Y0.X_weighted)
    weighted_MuDifference <- matrix(weighted_Mu1 - weighted_Mu0, ncol = 1)

    WSMDval <- sqrt(t(weighted_MuDifference) %*%
                      solve(var(Z1_NoME)) %*%
                      weighted_MuDifference)

    return(WSMD = WSMDval)
  }

  # Approximate WSMD Calculation Function----
  approxWSMDCal <- function(Weights) {
    # Weights <- Weights_BA

    approxMuDifference <- Sigma_eps %*% par_BA
    approxWSMD <- t(approxMuDifference) %*%
      solve(True_Cov_trt) %*%         # Use the true variance matrix
      approxMuDifference
    approxWSMD <- sqrt(approxWSMD)


    return(approxWSMD)
  }

  # Basic Information----
  # covariates
  ## Dimension of all covariates
  p_Z <- 2
  ## Dimension of error-prone covariates
  p_X <- 1
  ## Dimension of error-free covariates
  p_U <- p_Z - p_X
  ## expectation of error-prone data
  Mean_X <- 5
  ## expectation of error-free data
  Mean_U <- 10
  ## Expectation vector of covariates
  Mu_Z <- c(Mean_X, Mean_U)
  sd_X <- 1
  sd_U1 <- 1
  # outcome variable
  Sigma_Y <- 4
  rho_Y <- 0
  # Constant treatment effect: ATE = ATT
  ATE <- 10
  # True outcome model parameters
  Meancoef <- c(210, 27.4, 13.7)
  # True logistic model parameter
  True.theta <- c(0.5, -3, 1.5)

  # Simulation Information ----
  set.seed(102)
  vec_N <- 50000
  vec_rhoZ <- 0.3
  vec_Vareps <- seq(0, 0.5, by = 0.05)
  N_sim <- 1000

  # Result matrix
  mat_BA_rhoZ03_N10w <- matrix(ncol = 12,
                               nrow = length(vec_N) *
                                 length(vec_Vareps) *
                                 length(vec_rhoZ))
  colnames(mat_BA_rhoZ03_N10w) <- c("N", "rhoZ", "Var_eps",
                                    paste0("theta", 1:p_Z), "ATT",
                                    paste0("X", 1:p_X, "ASMD"), paste0("U", 1:p_U, "ASMD"),
                                    paste0("X", 1:p_X, "approxASMD"), paste0("U", 1:p_U, "approxASMD"),
                                    "SampleWSMD", "ApproxWSMD")

  # Use the "list" to save the result data of each simulation
  Coefdf_theta_BA_rhoZ03_N5kSim5 <- list()
  Biasdf_ATT_BA_rhoZ03_N5kSim5 <- list()
  df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5 <- list()
  df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5 <- list()
}

  # Progress bar
  pg <- 0
  # Counting row number for saving this time simulation: N x rhoZ x Var(eps)
  rowcount <- 0
  for (iiii in 1:N_sim) {
    # iiii: indicator of the iiii(th) time of simulation
    ## iiii <- 1

    for (iii in 1:length(vec_N)) {
      # iii: indicator of the iii(th) N
      ## iii <- 1
      ### N <- 100
      N <- vec_N[iii]

      for (ii in 1:vec_rhoZ) {
        # ii: indicator of the ii(th) rho_Z
        ## ii <- 1
        rhoZ <- vec_rhoZ[ii]

        for (i in 1:length(vec_Vareps)) {
          # i: indicator of the i(th) Var(eps)
          ## i <- 5
          Var_eps <- vec_Vareps[i]

          # Covariance matrix of covariates
          Sigma_Z <- matrix(c(sd_X^2, rhoZ * sd_X * sd_U1,
                              rhoZ * sd_X * sd_U1, sd_U1^2),
                            ncol = p_Z,
                            byrow = T)
          Sigma_eps <- diag(x = c(rep(Var_eps, p_X), rep(0, p_U)))

          # Covariates----
          ## ·Z----
          Z <- mvrnorm(n = N, mu = Mu_Z, Sigma = Sigma_Z)
          eps_X <- rnorm(n = N, mean = 0, sd = sqrt(Var_eps))
          colnames(Z) <- c(paste0("X", 1:p_X), paste0("U", 1:p_U))
          ### error-prone covariate
          X <- Z[, paste0("X", 1:p_X)]
          obsX <- X + eps_X
          ## ·Propensity score & # Treatment assignment----
          eta <- Z %*% True.theta[-1] + True.theta[1]
          pi <- exp(eta) / (1 + exp(eta))
          W <- sapply(X = pi, FUN = function(x) { rbinom(1, size = 1, x) })
          # sum(W)/N
          ## ·Outcome variable----
          # Potential outcomes
          ## Potential Outcome Y(0)
          Mean_i_0 <- Z %*% Meancoef[-1] + Meancoef[1]
          ## Potential Outcome Y(1) = Y(0) + constant treatment effect
          Mean_i_1 <- Mean_i_0 + ATE
          ## Potential Outcome Matrix: [Y(0), Y(1)]
          Mean_i_01 <- as.matrix(cbind(Mean_i_0, Mean_i_1))
          # Observed treatment effect
          Sigma_Y01 <- matrix(data = c(Sigma_Y, rho_Y * Sigma_Y,
                                       rho_Y * Sigma_Y, Sigma_Y),
                              ncol = 2, nrow = 2)
          eval <- eigen(Sigma_Y01, symmetric = TRUE)
          y_init <- matrix(stats::rnorm(N * 2, mean = 0, sd = 1), nrow = N, ncol = 2)
          y_tmp <- t(eval$vectors %*%
                       diag(sqrt(eval$values), nrow = 2) %*%
                       t(y_init))
          y_pot <- y_tmp + Mean_i_01
          Y <- (1 - W) * y_pot[, 1] + W * y_pot[, 2]
          # ·Data into frame ----
          Dat <- as.data.frame(cbind(Y, W, Z, obsX, y_pot))
          colnames(Dat) <- c("Y", "W", paste0("U", 1:p_Z), paste0("X", 1:p_X), paste0("Yi", c(0, 1)))
          ## Observed data
          W <- Dat$W
          Index_col <- which(W == 0)
          Index_trt <- which(W == 1)
          Y <- Dat$Y
          # Y0 <- Y[W == 0]
          Y0 <- Y[Index_col]
          # Y1 <- Y[W == 1]
          Y1 <- Y[Index_trt]

          ## Observed Covariates Z*
          Z <- Dat[, c(paste0("X", 1:p_X), paste0("U", p_X + 1:p_U))]
          Z0 <- as.matrix(Z[Index_col,])
          Z1 <- as.matrix(Z[Index_trt,])
          Z1_bar <- colMeans(Z1)
          Z <- as.matrix(Z)

          ## Observed Ture Covariate: Z without measurement error
          Z_NoME <- Dat[, paste0("U", 1:p_Z)]
          Z0_NoME <- as.matrix(Z_NoME[Index_col,])
          Z1_NoME <- as.matrix(Z_NoME[Index_trt,])
          Z1_bar_NoME <- colMeans(Z1_NoME)
          Z_NoME <- as.matrix(Z_NoME)

          # ZQY
          opt.out.zqy <- optim(par = rep(0, p_Z),
                               fn = objective.zqy,
                               gr = gradient.zqy,
                               method = "BFGS",
                               control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200))
          theta_star <- opt.out.zqy$par

          # ATT Estimation with theta*
          par_BA <- as.matrix(theta_star, ncol = 1)
          weight_BA <- exp(Z0_NoME %*% par_BA)
          weight_BA <- weight_BA / sum(weight_BA)
          ATT_BA <- mean(Y1) - t(weight_BA) %*% Y0
          ## Weights of theta* for calculating ASMD
          Weights_BA <- rep(1 / sum(W), N)
          Weights_BA[W == 0] <- weight_BA

          ASMD_BA <- ASMDCal(Weights = Weights_BA)
          approxASMD_BA <- approxASMDCal(Weights = Weights_BA)
          WSMD_BA <- WSMDCal(Weights = Weights_BA)
          approxWSMD_BA <- approxWSMDCal(Weights = Weights_BA)

          rowcount <- rowcount + 1
          mat_BA_rhoZ03_N10w[rowcount,] <- c(N, rhoZ, Var_eps,
                                             theta_star, ATT_BA,
                                             as.numeric(ASMD_BA),
                                             as.numeric(approxASMD_BA),
                                             WSMD_BA,
                                             approxWSMD_BA)

          pg <- pg + 1
          print(paste0("N = ", N, ", rhoZ = ", rhoZ, ", Var_eps = ", Var_eps))
          print(paste0("Progress = ", round(pg / (1 * length(vec_Vareps) * 1 * N_sim) * 100, digits = 2), "%"))
        }
      }
    }

    # Bias Data Frame Construction
    benchmark <- c(0, 0, 0,
                   0, 0, ATE,
                   0, 0, 0, 0,
                   0, 0)
    ### rhoZ = 0.3
    # by row: mat_BA_rhoZ03_N10w - benchmark
    Biasmat_BA_rhoZ03_N10w <- sweep(mat_BA_rhoZ03_N10w, MARGIN = 2, benchmark)

    Coefmat_BA <- Biasmat_BA_rhoZ03_N10w[, c("N", "Var_eps", "rhoZ",
                                             paste0("theta", 1:p_Z))]
    Biasmat_BA_ATT <- Biasmat_BA_rhoZ03_N10w[, c("N", "Var_eps", "rhoZ",
                                                 "ATT")]
    ## mat_ASMD and df_ASMD
    mat_BA_X1U1SampleASMD <- mat_BA_rhoZ03_N10w[, c("N", "Var_eps", "rhoZ",
                                                    paste0("X", 1:p_X, "ASMD"), paste0("U", 1:p_U, "ASMD"))]
    mat_BA_X1U1SampleASMD <- cbind(mat_BA_X1U1SampleASMD, "Sample")
    colnames(mat_BA_X1U1SampleASMD) <- c("N", "Var_eps", "rhoZ", "X1-SampleASMD", "U1-SampleASMD", "ASMDType")
    mat_BA_X1U1ApproxASMD <- mat_BA_rhoZ03_N10w[, c("N", "Var_eps", "rhoZ",
                                                    paste0("X", 1:p_X, "approxASMD"), paste0("U", 1:p_U, "approxASMD"))]
    mat_BA_X1U1ApproxASMD <- cbind(mat_BA_X1U1ApproxASMD, "Approx")
    colnames(mat_BA_X1U1ApproxASMD) <- c("N", "Var_eps", "rhoZ", "X1-ApproxASMD", "U1-ApproxASMD", "ASMDType")
    mat_BA_X1U1SampleApproxASMD <- rbind(mat_BA_X1U1SampleASMD, mat_BA_X1U1ApproxASMD)

    # 注意！melt的对象一定要是data.frame，否则无法melt成功
    Coefdf_theta_BA_rhoZ03_N5kSim5[[iiii]] <- melt(as.data.frame(Coefmat_BA),
                                                   id.vars = c("N", "Var_eps", "rhoZ"),
                                                   value.name = "Bias")
    Biasdf_ATT_BA_rhoZ03_N5kSim5[[iiii]] <- melt(as.data.frame(Biasmat_BA_ATT),
                                                 id.vars = c("N", "Var_eps", "rhoZ"),
                                                 value.name = "Bias")
    df_X1U1SampleASMD_BA_rhoZ03_N10w <- melt(as.data.frame(mat_BA_X1U1SampleASMD),
                                             id.vars = c("N", "Var_eps", "rhoZ", "ASMDType"),
                                             value.name = "ASMD")
    df_X1U1ApproxASMD_BA_rhoZ03_N10w <- melt(as.data.frame(mat_BA_X1U1ApproxASMD),
                                             id.vars = c("N", "Var_eps", "rhoZ", "ASMDType"),
                                             value.name = "ASMD")
    df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]] <- rbind(df_X1U1SampleASMD_BA_rhoZ03_N10w,
                                                               df_X1U1ApproxASMD_BA_rhoZ03_N10w)
    df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]]$ASMD <- as.numeric(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]]$ASMD)
    # class(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]]$ASMD)
    df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]]$Var_eps <- as.numeric(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]]$Var_eps)
    # class(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]]$Var_eps)
    df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]] <- cbind(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]],
                                                               rep(rep(c("X1", "U1"), each = nrow(mat_BA_X1U1SampleASMD)), 2))
    colnames(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[iiii]]) <- c("N", "Var_eps", "rhoZ", "ASMDType", "variable", "ASMD", "CovariateType")
    ## mat_WSMD and df_WSMD
    mat_BA_X1U1SampleApproxWSMD <- mat_BA_rhoZ03_N10w[, c("N", "Var_eps", "rhoZ",
                                                          "SampleWSMD", "ApproxWSMD")]
    # data.frame
    df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5[[iiii]] <- melt(as.data.frame(mat_BA_X1U1SampleApproxWSMD),
                                                        id.vars = c("N", "Var_eps", "rhoZ"),
                                                        value.name = "WSMD")
    # class(df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5$WSMD)

    # Result matrix
    mat_BA_rhoZ03_N10w <- matrix(ncol = 12,
                                 nrow = length(vec_N) *
                                   length(vec_Vareps) *
                                   length(vec_rhoZ))
    colnames(mat_BA_rhoZ03_N10w) <- c("N", "rhoZ", "Var_eps",
                                      paste0("theta", 1:p_Z), "ATT",
                                      paste0("X", 1:p_X, "ASMD"), paste0("U", 1:p_U, "ASMD"),
                                      paste0("X", 1:p_X, "approxASMD"), paste0("U", 1:p_U, "approxASMD"),
                                      "SampleWSMD", "ApproxWSMD")
    rowcount <- 0
  }


  # Data average
{
  variable <- Coefdf_theta_BA_rhoZ03_N5kSim5[[1]]$variable
  Coefdf_theta_BA_rhoZ03_N5kSim5 <- Reduce("+", Coefdf_theta_BA_rhoZ03_N5kSim5) / N_sim
  Coefdf_theta_BA_rhoZ03_N5kSim5$variable <- variable


  variable <- Biasdf_ATT_BA_rhoZ03_N5kSim5[[1]]$variable
  Biasdf_ATT_BA_rhoZ03_N5kSim5 <- Reduce("+", Biasdf_ATT_BA_rhoZ03_N5kSim5) / N_sim
  Biasdf_ATT_BA_rhoZ03_N5kSim5$variable <- variable

  variable <- df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[1]]
  meanASMD <- rep(0, length(variable$ASMD))
  for (i in 1:N_sim) {
    meanASMD <- meanASMD + df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5[[i]]$ASMD
  }
  meanASMD <- meanASMD / N_sim
  variable$ASMD <- meanASMD
  df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5 <- variable

  variable <- df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5[[1]]$variable
  df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5 <- Reduce("+", df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5) / N_sim
  df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5$variable <- variable
  # class(df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5)
}


  # Plot
{
  # Plot: N x rho
  ## Bias polt of coefficients----
  # summary(Coefdf_theta_BA_rhoZ03_N5kSim5$Bias)
  p_theta_rhoZ03_N5kSim5 <- ggplot(Coefdf_theta_BA_rhoZ03_N5kSim5) +
    geom_line(aes(x = Var_eps, y = Bias, color = variable)) +
    geom_point(aes(x = Var_eps, y = Bias, color = variable, shape = variable)) +
    # facet_grid(rows = vars(N), cols = vars(rhoZ)) + # 行按N进行分区，列按rhoZ分区
    labs(x = "Variance of Measurement Error", y = "Estimated Coefficients") +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
    scale_y_continuous(breaks = seq(-3.5, 2, by = 0.5)) +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.2) +
    guides(colour = guide_legend(title = NULL), shape = guide_legend(title = NULL))
  p_theta_rhoZ03_N5kSim5

  ## Absolute bias polt of ATT----
  # summary(abs(Biasdf_ATT_BA_rhoZ03_N5kSim5$Bias))
  maxval <- ceiling(summary(abs(Biasdf_ATT_BA_rhoZ03_N5kSim5$Bias))[6])
  pbias_ATT_rhoZ03_N5kSim5 <- ggplot(Biasdf_ATT_BA_rhoZ03_N5kSim5) +
    geom_line(aes(x = Var_eps, y = abs(Bias), color = variable)) +
    geom_point(aes(x = Var_eps, y = abs(Bias), color = variable, shape = variable),
               shape = 15, size = 2) +
    labs(x = "Variance of Measurement Error", y = "Absolute Bias of ATT") +
    theme_bw() +
    scale_y_continuous(breaks = seq(0, maxval, by = 1), limits = c(0, maxval)) +
    scale_x_continuous(breaks = seq(0, 0.6, by = 0.05)) +
    scale_color_manual(values = "#f1a340") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    guides(colour = guide_legend(title = NULL), shape = guide_legend(title = NULL))
  pbias_ATT_rhoZ03_N5kSim5


  ## ASMD plot：Contain both sample ASMD and approximate ASMD----
  # summary(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5$ASMD)
  p_X1U1SampleApproxASMD_rhoZ03_N5kSim5 <- ggplot(df_X1U1SampleApproxASMD_BA_rhoZ03_N5kSim5) +
    geom_line(aes(x = Var_eps, y = ASMD,
                  color = factor(x = ASMDType),
                  linetype = ASMDType,
                  group = variable),
              size = 0.8) +
    geom_point(aes(x = Var_eps, y = ASMD,
                   color = factor(x = ASMDType),
                   shape = CovariateType),
               size = 3.5) +
    scale_colour_manual(values = c("#7fbc41", "#c51b7d"), name = "ASMD Type") +
    scale_linetype(name = "ASMD Type") +
    scale_shape_manual(values = c(8, 18), name = "Covariate") +
    scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
    # 设置y-axis刻度线间隔
    scale_y_continuous(breaks = seq(0, 1.9, by = 0.1)) +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw() +
    labs(x = "Variance of Measurement Error", y = "Absolute Standardized Mean Difference (ASMD)") +
    theme(legend.position = "top")
  p_X1U1SampleApproxASMD_rhoZ03_N5kSim5
  override.shape <- c(".", ".")
  override.linetype <- c("solid", "dotted")
  p_X1U1SampleApproxASMD_rhoZ03_N5kSim5 <- p_X1U1SampleApproxASMD_rhoZ03_N5kSim5 +
    guides(colour = guide_legend(override.aes = list(shape = override.shape,
                                                     linetype = override.linetype)))
  p_X1U1SampleApproxASMD_rhoZ03_N5kSim5

  ## WSMD plot：Contain both sample WSMD and approximate WSMD----
  # summary(df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5$WSMD)
  df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5 <- df_X1U1SampleWSMD_BA_rhoZ03_N5kSim5
  colnames(df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5) <- c("N", "Var_eps", "rhoZ", "WSMDType", "WSMD")
  df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType <- as.character(df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType)
  RowIndex_Sample <- which(df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType == "SampleWSMD")
  df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType[RowIndex_Sample] <- "Sample"
  RowIndex_Approx <- which(df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType == "ApproxWSMD")
  df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType[RowIndex_Approx] <- c("Approx")
  df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType <- as.factor(df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5$WSMDType)

  p_X1U1SampleApproxWSMD_rhoZ03_N5kSim5 <- ggplot(df_X1U1SampleApproxWSMD_BA_rhoZ03_N5kSim5) +
    geom_line(aes(x = Var_eps, y = WSMD, color = WSMDType, linetype = WSMDType),
              size = 1) +
    geom_point(aes(x = Var_eps, y = WSMD, color = WSMDType, shape = WSMDType),
               size = 2) +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
    scale_y_continuous(breaks = seq(0, 1.56, by = 0.1)) +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw() +
    labs(x = "Variance of Measurement Error", y = "Mahalanobis Distance (MD)") +
    theme(legend.position = "top")
  p_X1U1SampleApproxWSMD_rhoZ03_N5kSim5

  Figure_BA <- ggarrange(p_theta_rhoZ03_N5kSim5, pbias_ATT_rhoZ03_N5kSim5,
                         p_X1U1SampleApproxASMD_rhoZ03_N5kSim5, p_X1U1SampleApproxWSMD_rhoZ03_N5kSim5,
                         nrow = 2, ncol = 2,
                         widths = c(8, 8, 8, 8), labels = c('(a)', '(b)', '(c)', '(d)'),
                         vjust = 1.2, hjust = 0,
                         common.legend = F, legend = 'top',
                         font.label = list(size = 15, face = 'plain'))
  Figure_BA

  filename_BA <- paste0("pBiasAnalysis_2Covariates_N", vec_N, "_Sim", N_sim, "_TrueVar.pdf")
  ggsave(filename = filename_BA, Figure_BA,
         width = 27, height = 27,
         units = "cm")
}

}