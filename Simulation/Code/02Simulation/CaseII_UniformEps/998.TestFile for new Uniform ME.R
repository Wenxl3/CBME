# -*- coding: utf-8 -*-
# ------------------------------------------
# Title       : 998.TestFile for new Uniform ME
# Remark      : The measurement error follows U(-1,1).
# ------------------------------------------

{
  rm(list = ls())
  getwd()
  # To suppress warnings in the global settings
  options(warn = -1)

  # Package Loading----
  library(MASS)
  library(reshape2)
  library(dplyr)
  library(gtools)
  library(ggplot2)
  library(stats)
  library(nleqslv)
  library(rootSolve)
  library(ggpubr)
  # Export the data of ASMD into the same Excel with different sheets
  library(openxlsx)

  # Importing Self-Defined Function----
  ## MSE for theta1, theta2, theta3, theta4, ATT
  MSE <- function(x) {
    # trueVal <- c(True.theta[-1], ATE)
    colMeans((t(t(x) - trueVal))^2)
  }

  # ASMD Calculation Function----
  ASMDCal <- function(Weights) {
    # Weights <- Weights_BA
    Weights_col <- Weights
    Weights <- W / sum(W)
    Weights[Weights == 0] <- Weights_col

    # Weighted Covariate (covariate without measurement error)
    WeightedCovariate <- Z * Weights

    # ASMD for each Covariate
    ## Numerator
    ASMD.Numerator_trt <- colSums(WeightedCovariate[Index_trt,])
    ASMD.Numerator_col <- colSums(WeightedCovariate[Index_col,])
    ## Denominator
    ASMD.Denominator1 <- apply(Z1, MARGIN = 2, sd)

    ##ASMD
    ASMD <- abs(ASMD.Numerator_trt - ASMD.Numerator_col) / ASMD.Denominator1

    return(data.frame(X1ASMD = ASMD[1], X2ASMD = ASMD[2],
                      U3ASMD = ASMD[3], U4ASMD = ASMD[4]))
  }

  ## Function_ObjectiveGradientHessian
  source(file = "./CEB/R/02Simulation/CaseII_UniformEps/999.ObjectiveGradient_UniformME.R")

  # Variable Information----
{
  ## Sample size
  N <- 2000
  ## numbers of replicate, i.e. N_Replicate
  r <- 2
  ## ·Covariates----
  ### X is error-prone data
  Mean_X <- c(4, 2)
  Mean_U <- c(3, 1)
  p_X <- length(Mean_X)
  p_U <- length(Mean_U)
  p_Z <- p_U + p_X
  ### Sigma_Z: Covariance matrix of covariates
  Sigma_Z <- diag(rep(1, 4))
  ### rhoZ: Correlation coefficient between error-prone covariates and error-free covariates
  # vec_rhoZ <- c(0, 3, 6) / 10
  vec_rhoZ <- 3 / 10
  ### Variance of measurement errors
  vec_Vareps <- seq(0, 0.5, by = 0.05)
  ## ·Propenstiy score model----
  ### True logistic model parameter
  True.theta <- c(3.5, -1, 0.5, -0.25, -0.1)
  ## ·Outcome variable----
  ### Diagonal element of the covariance matrix of outcome variables
  Sigma_Y <- 4
  ### True outcome model parameters
  Meancoef <- c(210, 27.4, 13.7, 13.7, 13.7)
  ### Correlation coefficient between Y(0) and Y(1)
  rho_Y <- 0
  ### True constant treatment effect
  ATE <- 10
}


  # List for replicated data----
{
  ## Uniform distributed measurement errors
  measure.error_Uniform <- list()
  X_Uniform <- list()
  Ze_Uniform <- list()
  Ze0_Uniform <- list()
  Ze1_Uniform <- list()
  Ze1_bar_Uniform <- list()
}

  # Simulation Setting----
{
  set.seed(seed = 102)
  ## N_Sim: Number of Monte Carlo experiments
  N_Sim <- 1000
}

  # Construct data frame for saving results----
{
{
  ### mat_SingleSim:
  #### 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
  #### nrow = N_Sim
  mat_SingleSim_eb_Uniform <- matrix(data = 0,
                                     ncol = 7, nrow = N_Sim)
  colnames(mat_SingleSim_eb_Uniform) <- c("rhoZ", "Var_eps",
                                          paste0("theta", 1:p_Z),
                                          "ATT")
  mat_SingleSim_ebme_Uniform <- mat_SingleSim_eb_Uniform
  mat_SingleSim_bc_Uniform <- mat_SingleSim_eb_Uniform
  mat_SingleSim_hcc_Uniform <- mat_SingleSim_eb_Uniform
  mat_SingleSim_hyj_Uniform <- mat_SingleSim_eb_Uniform

  ### mat_avgSim:
  #### 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
  #### nrow = length(vec_rhoZ) * length(vec_Vareps)
  mat_avgSim_eb_Uniform <- matrix(ncol = 7,
                                  nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(mat_avgSim_eb_Uniform) <- c("rhoZ", "Var_eps",
                                       paste0("theta", 1:p_Z),
                                       "ATT")
  mat_avgSim_ebme_Uniform <- mat_avgSim_eb_Uniform
  mat_avgSim_bc_Uniform <- mat_avgSim_eb_Uniform
  mat_avgSim_hcc_Uniform <- mat_avgSim_eb_Uniform
  mat_avgSim_hyj_Uniform <- mat_avgSim_eb_Uniform

  ### mat_BiasSdMSE:
  #### 17 cols: rhoZ, Var(eps), bias(theta1, theta2, theta3, theta4, ATT), sd(theta1, theta2, theta3, theta4, ATT), MSE(theta1, theta2, theta3, theta4, ATT)
  #### nrow = length(vec_rhoZ) * length(vec_Vareps)
  mat_BiasSdMSE_eb_Uniform <- matrix(ncol = 17,
                                     nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(mat_BiasSdMSE_eb_Uniform) <- c("rhoZ", "Var_eps",
                                          paste0("BIAS_", c(paste0("theta", 1:4), "ATT")),
                                          paste0("SD_", c(paste0("theta", 1:4), "ATT")),
                                          paste0("MSE_", c(paste0("theta", 1:4), "ATT")))
  mat_BiasSdMSE_ebme_Uniform <- mat_BiasSdMSE_eb_Uniform
  mat_BiasSdMSE_bc_Uniform <- mat_BiasSdMSE_eb_Uniform
  mat_BiasSdMSE_hcc_Uniform <- mat_BiasSdMSE_eb_Uniform
  mat_BiasSdMSE_hyj_Uniform <- mat_BiasSdMSE_eb_Uniform

  ## Iteration numbers and indicator pool
  ### EBME-BFGS
  iterbfgsind.r_Uniform <- rep(TRUE, N_Sim)
  ### BC - approximate bias corrector
  iterbcind.r_Uniform <- iterbfgsind.r_Uniform
  ### HCC
  iterhccind.r_Uniform <- iterbfgsind.r_Uniform
  ### HYJ
  iterhyjind.r_Uniform <- iterbfgsind.r_Uniform

  ## Rate of model convergence: Calculate the rate of convergence in each combination of rho_Z and var_eps
  ### 3 cols: rhoZ, Var(eps), average value of the rate of convergence
  ### nrow <- length(vec_rhoZ) * length(vec_Vareps)
  Rc_ebme_Uniform <- matrix(data = 0,
                            ncol = 3,
                            nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(Rc_ebme_Uniform) <- c("rhoZ", "Var_eps", "Rc")
  Rc_bc_Uniform <- Rc_ebme_Uniform
  Rc_hcc_Uniform <- Rc_ebme_Uniform
  Rc_hyj_Uniform <- Rc_ebme_Uniform

  # matrix markdown the ASMD
  mat_ASMD_CEB_Uniform <- NULL
  mat_ASMD_naiveEB_Uniform <- NULL
  mat_ASMD_BCEB_Uniform <- NULL
  mat_ASMD_HL_Uniform <- NULL
  mat_ASMD_HW_Uniform <- NULL
}

}

  # Simulation progress counting----
  pg <- 0
  Count_Row <- 0
  ## The ith combination of 1-RR and rhoZ
  ith_vareps_rhoZ <- 0
  mat_ASMD <- NULL
}

# Simulation----
for (rhoZ in vec_rhoZ) {
  # # Correlation coefficient between X and U
  # rhoZ <- 0.3

  # Add correlation to the covariance matrix of covariates
{
  Sigma_Z[1, 3] <- rhoZ
  Sigma_Z[3, 1] <- rhoZ
  Sigma_Z[2, 3] <- rhoZ
  Sigma_Z[3, 2] <- rhoZ
}

  for (var_eps in vec_Vareps) {
    # # Varaince of measurement error
    # var_eps <- 0.05

    # Covariance matrix of measurement errors containing only error-prone covariates.
    Sigma_e <- diag(rep(var_eps, p_X))
    ## Sigma: the covariance matrix of measurement errorss containing both error-prone and error-free covariates.
    Sigma <- diag(c(diag(Sigma_e), rep(0, p_U)))

    # mat_SingleSim_ASMD_CEB: markdown ASMD in each simulation
    ## ncol = dim(Z)
    mat_SingleSim_ASMD_CEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_naiveEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_BCEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HL <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HW <- matrix(0, nrow = N_Sim, ncol = p_Z)

    for (i in 1:N_Sim) {
      # i <- 1

      # Mark the rhoZ
    {
      ## Uniform distributed measurement errors
      mat_SingleSim_eb_Uniform[i, "rhoZ"] <- rhoZ
      mat_SingleSim_ebme_Uniform[i, "rhoZ"] <- rhoZ
      mat_SingleSim_bc_Uniform[i, "rhoZ"] <- rhoZ
      mat_SingleSim_hcc_Uniform[i, "rhoZ"] <- rhoZ
      mat_SingleSim_hyj_Uniform[i, "rhoZ"] <- rhoZ
    }

      # Mark the Var_eps
    {
      ## Uniform distributed measurement errors
      mat_SingleSim_eb_Uniform[i, "Var_eps"] <- var_eps
      mat_SingleSim_ebme_Uniform[i, "Var_eps"] <- var_eps
      mat_SingleSim_bc_Uniform[i, "Var_eps"] <- var_eps
      mat_SingleSim_hcc_Uniform[i, "Var_eps"] <- var_eps
      mat_SingleSim_hyj_Uniform[i, "Var_eps"] <- var_eps
    }

      # True covariates
    {
      Z <- mvrnorm(n = N, mu = c(Mean_X, Mean_U), Sigma = Sigma_Z)
      trueX <- Z[, 1:p_X]
      U <- Z[, (p_X + 1):p_Z]
    }

      # Construct replicated covariate dataset
    {
      # Uniform distributed measurement errors for replicated dataset
      ## Uniform(-1,1) measurement error

      # shape_1 indicates the lower bound of the uniform distribution
      shape_1 <- -1
      # shape_2 indicates the upper bound of the uniform distribution
      shape_2 <- 1

      # mean value of the uniform distribution
      mean_val_uniform <- (shape_1 + shape_2) / 2
      # variance value of the uniform distribution
      var_val_uniform <- (shape_2 - shape_1)^2 / 12

      for (k in 1:r) {
        if (var_eps == 0) {
          measure.error_Uniform[[k]] <- rep(0, p_X * N)
        }else {
          # measurement error follows U[shape_1, shape_2]
          measure.error_Uniform[[k]] <- runif(n = p_X * N, min = shape_1, max = shape_2)

          # Adjust the expectation of measurement error to 0
          measure.error_Uniform[[k]] <- measure.error_Uniform[[k]] - mean_val_uniform
          # Adjust the variance of measurement error to prespecified variance
          measure.error_Uniform[[k]] <- measure.error_Uniform[[k]] * sqrt(var_eps) / sqrt(var_val_uniform)
        }

        measure.error_Uniform[[k]] <- matrix(data = measure.error_Uniform[[k]], nrow = N, ncol = p_X)

        X_Uniform[[k]] <- trueX + measure.error_Uniform[[k]]
        Ze_Uniform[[k]] <- data.frame(cbind(X_Uniform[[k]], U))
        colnames(Ze_Uniform[[k]]) <- c("X1", "X2", "U3", "U4")
      }
    }

      # Propensity score
    {
      eta <- Z %*% True.theta[-1] + True.theta[1]
      pi <- exp(eta) / (1 + exp(eta))
      # Treatment assignment
      W <- sapply(X = pi, FUN = function(x) { rbinom(1, size = 1, x) })
      Index_col <- which(W == 0)
      Index_trt <- which(W == 1)
    }

      # Outcome variable
    {
      # treatment effect
      Mean_i_0 <- Z %*% Meancoef[-1] + Meancoef[1]
      Mean_i_1 <- Mean_i_0 + ATE
      Mean_i_01 <- as.matrix(cbind(Mean_i_0, Mean_i_1))
      # Potential outcomes
      Sigma_Y01 <- matrix(data = c(Sigma_Y, rho_Y * Sigma_Y, rho_Y * Sigma_Y, Sigma_Y),
                          ncol = 2, nrow = 2)
      eval <- eigen(Sigma_Y01, symmetric = TRUE)
      y_init <- matrix(stats::rnorm(N * 2, 0, 1), nrow = N, ncol = 2)
      y_tmp <- t(eval$vectors %*%
                   diag(sqrt(eval$values), nrow = 2) %*%
                   t(y_init))
      y_pot <- y_tmp + Mean_i_01
      Y <- (1 - W) * y_pot[, 1] + W * y_pot[, 2]
      Y0 <- Y[Index_col]
      Y1 <- Y[Index_trt]
    }

      # True covariates Z
    {
      Z0 <- as.matrix(Z[Index_col,])
      Z1 <- as.matrix(Z[Index_trt,])
      Z1_bar <- colMeans(Z1)
    }

      # Mean value of observed Z1 (with measurement errors) in replicated dataset
    {
      Ze1_bar_al_Uniform <- rep(0, p_Z)
      for (k in 1:r) {
        Ze0_Uniform[[k]] <- as.matrix(Ze_Uniform[[k]][Index_col,])
        Ze1_Uniform[[k]] <- as.matrix(Ze_Uniform[[k]][Index_trt,])
        Ze1_bar_Uniform[[k]] <- colMeans(Ze1_Uniform[[k]])
        Ze1_bar_al_Uniform <- Ze1_bar_al_Uniform + Ze1_bar_Uniform[[k]]
        Ze_Uniform[[k]] <- as.matrix(Ze_Uniform[[k]])
      }
      Ze1_bar_al_Uniform <- Ze1_bar_al_Uniform / r
    }

      # Estimate the variance-covariance matrix $\Sigma$ of measurement error
    {
      meanZe_al_Uniform <- Reduce("+", Ze_Uniform) / r
      erhat_Uniform <- matrix(rep(0, ncol(Z)^2), ncol = ncol(Z), nrow = ncol(Z))
      for (kk in 1:r) {
        erhat_Uniform <- erhat_Uniform +
          t(Ze_Uniform[[kk]] - meanZe_al_Uniform) %*% (Ze_Uniform[[kk]] - meanZe_al_Uniform)
      }
      Sigmahat_al_Uniform <- erhat_Uniform / ((r - 1) * N)
      # round(Sigmahat_al_Uniform, 3)
    }

      # ZQY-EB: native Entropy Balancing
    {
      opt.out.zqy.r_Uniform <- optim(par = rep(0, (ncol(Z))),
                                     fn = objective.zqy.r_Uniform,
                                     gr = gradient.zqy.r_Uniform,
                                     method = "BFGS",
                                     control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200))
      theta0_Uniform <- opt.out.zqy.r_Uniform$par
      # round(theta0_Uniform, 3)
    }

      # EBME - BFGS
    {
      output.CEB.solveObj_Uniform <- tryCatch(optim(par = theta0_Uniform,
                                                    fn = objective.ebme.r_Uniform,
                                                    gr = gradient.ebme.r_Uniform,
                                                    method = "BFGS",
                                                    control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200)),
                                              error = function(e) {
                                                print("Fail to solve objective of CEB, We directly use naive EB's results for subsequent comparison.")
                                                # 此时无法求解CEB，直接将solve CEB Objective的结果赋为FALSE
                                                return(FALSE)
                                              })
      # 如果optim无法求解，output.CEB.solveObj_Uniform会是一个logical变量；否则，如果optim能求解，output.CEB.solveObj_Uniform是一个list变量
      if (class(output.CEB.solveObj_Uniform) != "list") {
        # 此时无法求解CEB，将naive EB的结果赋给CEB
        output.CEB.solveObj_Uniform <- opt.out.zqy.r_Uniform
      }
      if (is.na(output.CEB.solveObj_Uniform$value) | is.nan(output.CEB.solveObj_Uniform$value)) {
        # 此时无法求解CEB，将naive EB的结果赋给CEB
        output.CEB.solveObj_Uniform <- opt.out.zqy.r_Uniform
      }
      theta_CEB_solveObj <- output.CEB.solveObj_Uniform$par
      ObjValue_CEB_solveObj <- output.CEB.solveObj_Uniform$value

      # theta_CEB - solve the gradient funtion
      output.CEB.solveGradient_Uniform <- tryCatch(multiroot(f = gradient.ebme.r_Uniform,
                                                             start = theta0_Uniform,
                                                             rtol = .Machine$double.eps, maxiter = 200),
                                                   error = function(e) {
                                                     print("Fail to solve gradient of CEB, We directly use naive EB's results for subsequent comparison.")
                                                     # 此时无法求解CEB，直接将solve CEB Objective的结果赋为FALSE
                                                     return(FALSE)
                                                   })
      # 如果multiroot无法求解，output.CEB.solveGradient会是一个logical变量；否则，如果multiroot能求解，output.CEB.solveGradient是一个list变量
      if (class(output.CEB.solveGradient_Uniform) != "list") {
        # 此时求解CEB gradient 不收敛，直接赋值9999
        ObjValue_CEB_solveGradient <- 9999
      }else {
        theta_CEB_solveGradient <- output.CEB.solveGradient_Uniform$root
        ObjValue_CEB_solveGradient <- objective.ebme.r_Uniform(theta = theta_CEB_solveGradient)
        if (is.nan(ObjValue_CEB_solveGradient) | is.na(ObjValue_CEB_solveGradient)) { ObjValue_CEB_solveGradient <- 9999 }
      }

      if (ObjValue_CEB_solveObj <= ObjValue_CEB_solveGradient) {
        theta_CEB_Uniform <- theta_CEB_solveObj
        iterbfgsind.r_Uniform[i] <- (max(abs(gradient.ebme.r_Uniform(theta_CEB_solveObj))) < 1e-04)
        ## double check
        if (is.na(iterbfgsind.r_Uniform[i]) | is.nan(iterbfgsind.r_Uniform[i])) {
          iterbfgsind.r_Uniform[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }else {
        theta_CEB_Uniform <- theta_CEB_solveGradient
        iterbfgsind.r_Uniform[i] <- (max(abs(output.CEB.solveGradient_Uniform$f.root)) < sqrt(.Machine$double.eps)) & (output.CEB.solveGradient_Uniform$estim.precis != "NaN")
        ## double check
        if (is.na(iterbfgsind.r_Uniform[i]) | is.nan(iterbfgsind.r_Uniform[i])) {
          iterbfgsind.r_Uniform[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }
    }

      # BCEB
    {
      theta_bc_Uniform <- drop(ginv(hessian.zqy.r_Uniform(theta0_Uniform) - Sigma) %*%
                                 hessian.zqy.r_Uniform(theta0_Uniform) %*%
                                 theta0_Uniform)

      ### 判断BC是否无法求解：（默认【收敛】）存在以下情况任意其一都被判为【无解/不收敛】
      iterbcind.r_Uniform[i] <- !(NaN %in% theta_bc_Uniform | (sum(is.infinite(theta_bc_Uniform)) > 0))
    }

      # HCC(2012)
    {
      opt.out.hcc.r_Uniform <- multiroot(f = gradient.hcc.r_Uniform, start = rep(0, p_Z),
                                         rtol = .Machine$double.eps, maxiter = 200)
      iterhccind.r_Uniform[i] <- (max(gradient.hcc.r_Uniform(opt.out.hcc.r_Uniform$root)) < sqrt(.Machine$double.eps)) & (opt.out.hcc.r_Uniform$estim.precis != "NaN")
    }

      # HYJ(2012)
    {
      opt.out.hyj.r_Uniform <- multiroot(f = gradient.hyj.r_Uniform, start = rep(0, p_Z),
                                         rtol = .Machine$double.eps, maxiter = 200)
      iterhyjind.r_Uniform[i] <- (max(gradient.hyj.r_Uniform(opt.out.hyj.r_Uniform$root)) < sqrt(.Machine$double.eps)) & (opt.out.hyj.r_Uniform$estim.precis != "NaN")
      ### double check
      if (is.na(iterhyjind.r_Uniform[i])) {
        iterhyjind.r_Uniform[i] <- FALSE
        print("NA clear!")
      }
    }

      # ATT estimation with \hat{\theta}
    {
      # ZQY-EB: naive Entropy Balancing
    {
      # par_zqy_Uniform <- as.matrix(theta0_Uniform, ncol = 1)
      # weight_zqy_Uniform <- exp(Ze0_Uniform[[1]] %*% par_zqy_Uniform)
      weight_zqy_Uniform <- exp(Ze0_Uniform[[1]] %*% theta0_Uniform)
      weight_zqy_Uniform <- weight_zqy_Uniform / sum(weight_zqy_Uniform)
      ATT_zqy_Uniform <- mean(Y1) - t(weight_zqy_Uniform) %*% Y0
      mat_SingleSim_eb_Uniform[i, "ATT"] <- mat_SingleSim_eb_Uniform[i, "ATT"] + ATT_zqy_Uniform
    }

      # CEB / EBME - bfgs C2: the proposed method
    {
      # par_bfgs_Uniform <- as.matrix(theta_CEB_Uniform, ncol = 1)
      # weight_bfgs_Uniform <- exp(Ze0_Uniform[[1]] %*% par_bfgs_Uniform)
      weight_bfgs_Uniform <- exp(Ze0_Uniform[[1]] %*% theta_CEB_Uniform)
      weight_bfgs_Uniform <- weight_bfgs_Uniform / sum(weight_bfgs_Uniform)
      ATT_bfgs_Uniform <- mean(Y1) - t(weight_bfgs_Uniform) %*% Y0
      mat_SingleSim_ebme_Uniform[i, "ATT"] <- mat_SingleSim_ebme_Uniform[i, "ATT"] + ATT_bfgs_Uniform
    }

      # BCEB
    {
      # par_bc_Uniform <- as.matrix(theta_bc_Uniform, ncol = 1)
      # weight_bc_Uniform <- exp(Ze0_Uniform[[1]] %*% par_bc_Uniform)
      weight_bc_Uniform <- exp(Ze0_Uniform[[1]] %*% theta_bc_Uniform)
      weight_bc_Uniform <- weight_bc_Uniform / sum(weight_bc_Uniform)
      ATT_bc_Uniform <- mean(Y1) - t(weight_bc_Uniform) %*% Y0
      mat_SingleSim_bc_Uniform[i, "ATT"] <- mat_SingleSim_bc_Uniform[i, "ATT"] + ATT_bc_Uniform
    }

      # HCC
    {
      par_hcc_Uniform <- as.matrix(opt.out.hcc.r_Uniform$root, ncol = 1)
      weight_hcc_Uniform <- matrix(rep(0, nrow(weight_bc_Uniform)))
      for (k in 1:r) {
        weight_hcc_Uniform <- weight_hcc_Uniform + exp(Ze0_Uniform[[k]] %*% par_hcc_Uniform)
      }
      weight_hcc_Uniform <- weight_hcc_Uniform / sum(weight_hcc_Uniform)
      ATT_hcc_Uniform <- mean(Y1) - t(weight_hcc_Uniform) %*% Y0
      mat_SingleSim_hcc_Uniform[i, "ATT"] <- mat_SingleSim_hcc_Uniform[i, "ATT"] + ATT_hcc_Uniform
    }

      # HYJ
    {
      par_hyj_Uniform <- as.matrix(opt.out.hyj.r_Uniform$root, ncol = 1)
      weight_hyj_Uniform <- matrix(rep(0, nrow(weight_bc_Uniform)))
      for (k in 1:r) {
        weight_hyj_Uniform <- weight_hyj_Uniform + exp(Ze0_Uniform[[k]] %*% par_hyj_Uniform)
      }
      weight_hyj_Uniform <- weight_hyj_Uniform / sum(weight_hyj_Uniform)
      ATT_hyj_Uniform <- mean(Y1) - t(weight_hyj_Uniform) %*% Y0
      ### double check
      if (is.na(ATT_hyj_Uniform)) {
        iterhyjind.r_Uniform[i] <- FALSE
        print("ATT of HW - NA clear!")
      }
      mat_SingleSim_hyj_Uniform[i, "ATT"] <- mat_SingleSim_hyj_Uniform[i, "ATT"] + ATT_hyj_Uniform
    }

    }

      # theta coefficient
    {
      mat_SingleSim_eb_Uniform[i, paste0("theta", 1:p_Z)] <- theta0_Uniform
      mat_SingleSim_ebme_Uniform[i, paste0("theta", 1:p_Z)] <- theta_CEB_Uniform
      mat_SingleSim_bc_Uniform[i, paste0("theta", 1:p_Z)] <- theta_bc_Uniform
      mat_SingleSim_hcc_Uniform[i, paste0("theta", 1:p_Z)] <- opt.out.hcc.r_Uniform$root
      mat_SingleSim_hyj_Uniform[i, paste0("theta", 1:p_Z)] <- opt.out.hyj.r_Uniform$root
    }

      # ASMD calculation ----
    {
      # mat_SingleSim_ASMD_CEB
      mat_SingleSim_ASMD_naiveEB[i,] <- unlist(ASMDCal(weight_zqy_Uniform))
      mat_SingleSim_ASMD_CEB[i,] <- unlist(ASMDCal(weight_bfgs_Uniform))
      mat_SingleSim_ASMD_BCEB[i,] <- unlist(ASMDCal(weight_bc_Uniform))
      mat_SingleSim_ASMD_HL[i,] <- unlist(ASMDCal(weight_hcc_Uniform))
      mat_SingleSim_ASMD_HW[i,] <- unlist(ASMDCal(weight_hyj_Uniform))
    }

      # Output program progress
      pg <- pg + 1
      print(paste0("rhoZ = ", rhoZ, ", var_eps = ", var_eps, ", Sim.i = ", i))
      print(paste0("Progress = ",
                   round(pg / (length(vec_Vareps) *
                     length(vec_rhoZ) *
                     N_Sim) * 100,
                         digits = 2), "%"))
    }

    ith_vareps_rhoZ <- ith_vareps_rhoZ + 1
    print(paste0("The ", ith_vareps_rhoZ, "th row of mat_avgSim value assignment -- start"))
    # Mark down the rate of convergence
  {
    ## CEB - EBME
    Rc_ebme_Uniform[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterbfgsind.r_Uniform) / N_Sim)
    ## BC
    Rc_bc_Uniform[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterbcind.r_Uniform) / N_Sim)
    ## HCC
    Rc_hcc_Uniform[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterhccind.r_Uniform) / N_Sim)
    ## HYJ
    Rc_hyj_Uniform[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterhyjind.r_Uniform) / N_Sim)
  }

    # mat_avgSim赋值 - Convergent results
  {
    mat_avgSim_eb_Uniform[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_eb_Uniform)
    mat_avgSim_ebme_Uniform[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_ebme_Uniform[iterbfgsind.r_Uniform,])
    mat_avgSim_bc_Uniform[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_bc_Uniform[iterbcind.r_Uniform,])
    mat_avgSim_hcc_Uniform[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_hcc_Uniform[iterhccind.r_Uniform,])
    mat_avgSim_hyj_Uniform[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_hyj_Uniform[iterhyjind.r_Uniform,])
  }

    # mat_BiasSDMSE赋值 - Convergent results
    trueVal <- c(True.theta[-1], ATE)
  {
    ## rhoZ赋值
    mat_BiasSdMSE_eb_Uniform[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_ebme_Uniform[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_bc_Uniform[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_hcc_Uniform[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_hyj_Uniform[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    ## Var_eps
    mat_BiasSdMSE_eb_Uniform[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_ebme_Uniform[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_bc_Uniform[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_hcc_Uniform[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_hyj_Uniform[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    ## bias
    ColnamesinNeed_avgmat <- c(paste0("theta", 1:p_Z), "ATT")
    ColnamesinNeed_Bias <- paste0("BIAS_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Uniform[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_eb_Uniform[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_ebme_Uniform[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_ebme_Uniform[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_bc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_bc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_hcc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_hcc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_hyj_Uniform[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_hyj_Uniform[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    ## sd
    ColnamesinNeed_SD <- paste0("SD_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Uniform[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_eb_Uniform[, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_ebme_Uniform[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_ebme_Uniform[iterbfgsind.r_Uniform, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_bc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_bc_Uniform[iterbcind.r_Uniform, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_hcc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_hcc_Uniform[iterhccind.r_Uniform, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_hyj_Uniform[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_hyj_Uniform[iterhyjind.r_Uniform, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    ## MSE
    ColnamesinNeed_MSE <- paste0("MSE_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Uniform[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_eb_Uniform[, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_ebme_Uniform[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_ebme_Uniform[iterbfgsind.r_Uniform, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_bc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_bc_Uniform[iterbcind.r_Uniform, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_hcc_Uniform[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_hcc_Uniform[iterhccind.r_Uniform, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_hyj_Uniform[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_hyj_Uniform[iterhyjind.r_Uniform, ColnamesinNeed_avgmat])
  }

    # ASMD
  {
    # mat_ASMD_CEB_Uniform
    mat_ASMD_CEB_Uniform <- rbind(mat_ASMD_CEB_Uniform,
                                  colMeans(mat_SingleSim_ASMD_CEB, na.rm = TRUE))
    mat_ASMD_naiveEB_Uniform <- rbind(mat_ASMD_naiveEB_Uniform,
                                      colMeans(mat_SingleSim_ASMD_naiveEB, na.rm = TRUE))
    mat_ASMD_BCEB_Uniform <- rbind(mat_ASMD_BCEB_Uniform,
                                   colMeans(mat_SingleSim_ASMD_BCEB, na.rm = TRUE))
    mat_ASMD_HL_Uniform <- rbind(mat_ASMD_HL_Uniform,
                                 colMeans(mat_SingleSim_ASMD_HL, na.rm = TRUE))
    mat_ASMD_HW_Uniform <- rbind(mat_ASMD_HW_Uniform,
                                 colMeans(mat_SingleSim_ASMD_HW, na.rm = TRUE))
    # round(mat_ASMD_CEB_Uniform, 3)
  }

    # 清空记录上一个var_eps × rhoZ组合的500次simulation结果的数据集
    Count_Row <- 0
  {
  {
    # mat_SingleSim: 记录var_eps × rhoZ的某一个组合下的N_Sim次simulation的结果
    ## 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
    ## nrow = N_Sim
    mat_SingleSim_eb_Uniform <- matrix(data = 0,
                                       ncol = 7, nrow = N_Sim)
    colnames(mat_SingleSim_eb_Uniform) <- c("rhoZ", "Var_eps",
                                            paste0("theta", 1:p_Z), "ATT")
    mat_SingleSim_ebme_Uniform <- mat_SingleSim_eb_Uniform
    mat_SingleSim_bc_Uniform <- mat_SingleSim_eb_Uniform
    mat_SingleSim_hcc_Uniform <- mat_SingleSim_eb_Uniform
    mat_SingleSim_hyj_Uniform <- mat_SingleSim_eb_Uniform

    # 清空记录上一个var_eps × rhoZ组合下500次simulation的收敛结果标记
    ## EBME-BFGS
    iterbfgsind.r_Uniform <- rep(TRUE, N_Sim)
    ## BC - approximate bias corrector
    iterbcind.r_Uniform <- iterbfgsind.r_Uniform
    ## HCC
    iterhccind.r_Uniform <- iterbfgsind.r_Uniform
    ## HYJ
    iterhyjind.r_Uniform <- iterbfgsind.r_Uniform

    # 矩阵清空为下一轮simulation做准备
    mat_SingleSim_ASMD_naiveEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_CEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_BCEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HL <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HW <- matrix(0, nrow = N_Sim, ncol = p_Z)
  }
  }
  }
}


# Save .Rdata
save.image(file = "CompareAllMethods_UniformEps_Original.RData")

{
  # Table output----
{
  Rc_ebme_Uniform <- cbind(Rc_ebme_Uniform, Method = "CEB")
  Rc_bc_Uniform <- cbind(Rc_bc_Uniform, Method = "BCEB")
  Rc_hcc_Uniform <- cbind(Rc_hcc_Uniform, Method = "HL")
  Rc_hyj_Uniform <- cbind(Rc_hyj_Uniform, Method = "HW")
  mat_Rc_Uniform <- rbind(Rc_ebme_Uniform,
                          Rc_bc_Uniform,
                          Rc_hcc_Uniform,
                          Rc_hyj_Uniform)
  csvName_mat_Rc_Uniform <- paste0("mat_UniformEps_ConvergenceRate_N", N,
                                   "_Sim", N_Sim, ".csv")
  write.table(mat_Rc_Uniform, file = csvName_mat_Rc_Uniform,
              sep = ",", row.names = F)

  ## mat_avgSim
  mat_avgSim_eb_Uniform <- cbind(mat_avgSim_eb_Uniform, Method = "Naive")
  mat_avgSim_ebme_Uniform <- cbind(mat_avgSim_ebme_Uniform, Method = "CEB")
  mat_avgSim_bc_Uniform <- cbind(mat_avgSim_bc_Uniform, Method = "BCEB")
  mat_avgSim_hcc_Uniform <- cbind(mat_avgSim_hcc_Uniform, Method = "HL")
  mat_avgSim_hyj_Uniform <- cbind(mat_avgSim_hyj_Uniform, Method = "HW")
  mat_avgSim_Uniform <- rbind(mat_avgSim_eb_Uniform,
                              mat_avgSim_ebme_Uniform,
                              mat_avgSim_bc_Uniform,
                              mat_avgSim_hcc_Uniform,
                              mat_avgSim_hyj_Uniform)
  csvName_mat_avgSim <- paste0("mat_UniformEps_avgSim_theta_ATT_N", N, "_Sim", N_Sim,
                               ".csv")
  write.table(mat_avgSim_Uniform, file = csvName_mat_avgSim,
              sep = ",", row.names = F)

  ## mat_BiasSdMSE
  mat_BiasSdMSE_eb_Uniform <- cbind(mat_BiasSdMSE_eb_Uniform, Method = "Naive")
  mat_BiasSdMSE_ebme_Uniform <- cbind(mat_BiasSdMSE_ebme_Uniform, Method = "CEB")
  mat_BiasSdMSE_bc_Uniform <- cbind(mat_BiasSdMSE_bc_Uniform, Method = "BCEB")
  mat_BiasSdMSE_hcc_Uniform <- cbind(mat_BiasSdMSE_hcc_Uniform, Method = "HL")
  mat_BiasSdMSE_hyj_Uniform <- cbind(mat_BiasSdMSE_hyj_Uniform, Method = "HW")
  mat_BiasSdMSE_Uniform <- rbind(mat_BiasSdMSE_eb_Uniform,
                                 mat_BiasSdMSE_ebme_Uniform,
                                 mat_BiasSdMSE_bc_Uniform,
                                 mat_BiasSdMSE_hcc_Uniform,
                                 mat_BiasSdMSE_hyj_Uniform)
  csvName_mat_BiasSdMSE <- paste0("mat_UniformEps_BiasSdMSE_theta_ATT_N", N, "_Sim", N_Sim,
                                  ".csv")
  write.table(mat_BiasSdMSE_Uniform, file = csvName_mat_BiasSdMSE,
              sep = ",", row.names = F)
}

  # Plots----
{
  # Plot of convergenc rate of EBME, BC, HCC, HYJ
  Rcdf_Uniform <- as.data.frame(mat_Rc_Uniform)
  # Rcdf_Uniform$Method <- factor(Rcdf_Uniform$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  Rcdf_Uniform$Method <- factor(Rcdf_Uniform$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  p_Rc_UniformEps <- ggplot(Rcdf_Uniform) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = Rc, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = Rc, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Rate of Convergence", title = "Measurement Error with a Uniform Distribution") +
    theme_bw()
  p_Rc_UniformEps
  ggsave(filename = "p_UniformEps_Rc.pdf", p_Rc_UniformEps,
         width = 14, height = 13,
         units = "cm")

  ## Seperate bias, sd, MSE dataframe
  mat_bias_Uniform <- mat_BiasSdMSE_Uniform[, c("rhoZ", "Var_eps", paste0("BIAS_", c(paste0("theta", 1:4), "ATT")), "Method")]
  mat_sd_Uniform <- mat_BiasSdMSE_Uniform[, c("rhoZ", "Var_eps", paste0("SD_", c(paste0("theta", 1:4), "ATT")), "Method")]
  mat_MSE_Uniform <- mat_BiasSdMSE_Uniform[, c("rhoZ", "Var_eps", paste0("MSE_", c(paste0("theta", 1:4), "ATT")), "Method")]
  ## bias
  biasdf_Uniform <- melt(as.data.frame(mat_bias_Uniform),
                         id.vars = c("rhoZ", "Var_eps", "Method"))
  biasdf_Uniform$value <- as.numeric(biasdf_Uniform$value)
  which(is.na(biasdf_Uniform$value))
  # biasdf_Uniform$Method <- factor(biasdf_Uniform$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  biasdf_Uniform$Method <- factor(biasdf_Uniform$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### bias_theta
  biasdf_theta_Uniform <- biasdf_Uniform[biasdf_Uniform$variable != "BIAS_ATT",]
  p_AllMethods_bias_theta_UniformEps <- ggplot(biasdf_theta_Uniform) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) + # 行按variable进行分区，列按rhoZ分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Bias", title = "True Value: theta1 = -1, theta2 = 0.5, theta2 = -0.25, theta4 = -0.1") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    # ylim(-0.8, 0.6) +
    # scale_y_continuous(breaks = seq(-0.8, 0.6, by = 0.1)) +
    theme_bw()
  p_AllMethods_bias_theta_UniformEps
  ggsave(filename = "p_UniformEps_AllMethods_bias_theta.pdf", p_AllMethods_bias_theta_UniformEps,
         width = 26, height = 26,
         units = "cm")
  ### bias_ATT
  biasdf_ATT_Uniform <- biasdf_Uniform[biasdf_Uniform$variable == "BIAS_ATT",]
  p_AllMethods_bias_ATT_UniformEps <- ggplot(biasdf_ATT_Uniform) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    # labs(x = "Variance of Measurement Error", y = "Bias", title = "Bias of ATT with Uniform Measurement Errors") +
    labs(x = "Variance of Measurement Error", y = "Bias") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    # ylim(-9, 2) +
    # scale_y_continuous(breaks = seq(-9, 2, by = 0.5)) +
    theme_bw()
  p_AllMethods_bias_ATT_UniformEps
  ggsave(filename = "p_UniformEps_AllMethods_bias_ATT.pdf", p_AllMethods_bias_ATT_UniformEps,
         width = 26, height = 10,
         units = "cm")
  ## sd
  sddf_Uniform <- melt(as.data.frame(mat_sd_Uniform),
                       id.vars = c("rhoZ", "Var_eps", "Method"))
  sddf_Uniform$value <- as.numeric(sddf_Uniform$value)
  which(is.na(sddf_Uniform$value))
  # sddf_Uniform$Method <- factor(sddf_Uniform$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  sddf_Uniform$Method <- factor(sddf_Uniform$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### sd_theta
  sddf_theta_Uniform <- sddf_Uniform[sddf_Uniform$variable != "SD_ATT",]
  p_AllMethods_sd_theta_UniformEps <- ggplot(sddf_theta_Uniform) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) + # 行按variable进行分区，列按rhoZ分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Standard Deviation", title = "SD of Coefficients with Normal Distributed Measurement Error") +
    theme_bw()
  p_AllMethods_sd_theta_UniformEps
  ggsave(filename = "p_UniformEps_AllMethods_sd_theta.pdf", p_AllMethods_sd_theta_UniformEps,
         width = 26, height = 26,
         units = "cm")
  ### sd_ATT
  sddf_ATT_Uniform <- sddf_Uniform[sddf_Uniform$variable == "SD_ATT",]
  p_AllMethods_sd_ATT_UniformEps <- ggplot(sddf_ATT_Uniform) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    # scale_y_continuous(limits = c(0, 2.25), breaks = seq(from = 0, to = 2.25, by = 0.5)) +
    scale_y_continuous(limits = c(0, max(sddf_ATT_Uniform$value)), breaks = seq(from = 0, to = ceiling(max(sddf_ATT_Uniform$value)), by = 0.5)) +
    # labs(x = "Variance of Measurement Error", y = "Standard Deviation", title = "SD of ATT with Normal Measurement Errors") +
    labs(x = "Variance of Measurement Error", y = "Standard Deviation") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_sd_ATT_UniformEps
  ggsave(filename = "p_UniformEps_AllMethods_sd_ATT.pdf", p_AllMethods_sd_ATT_UniformEps,
         width = 26, height = 10,
         units = "cm")
  ## MSE
  MSEdf_Uniform <- melt(as.data.frame(mat_MSE_Uniform),
                        id.vars = c("rhoZ", "Var_eps", "Method"))
  MSEdf_Uniform$value <- as.numeric(MSEdf_Uniform$value)
  which(is.na(MSEdf_Uniform$value))
  # MSEdf_Uniform$Method <- factor(MSEdf_Uniform$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  MSEdf_Uniform$Method <- factor(MSEdf_Uniform$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### MSE_theta
  MSEdf_theta_Uniform <- MSEdf_Uniform[MSEdf_Uniform$variable != "MSE_ATT",]
  p_AllMethods_MSE_theta_UniformEps <- ggplot(MSEdf_theta_Uniform) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Mean Squared Error", title = "MSE of Coefficients with Uniform Distributed Measurement Error") +
    theme_bw()
  p_AllMethods_MSE_theta_UniformEps
  ggsave(filename = "p_UniformEps_AllMethods_MSE_theta.pdf", p_AllMethods_MSE_theta_UniformEps,
         width = 26, height = 26,
         units = "cm")
  ### MSE_ATT
  MSEdf_ATT_Uniform <- MSEdf_Uniform[MSEdf_Uniform$variable == "MSE_ATT",]
  p_AllMethods_MSE_ATT_UniformEps <- ggplot(MSEdf_ATT_Uniform) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    # labs(x = "Variance of Measurement Error", y = "Mean Squared Error", title = "MSE of ATT with Uniform Measurement Errors") +
    labs(x = "Variance of Measurement Error", y = "Mean Squared Error") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_MSE_ATT_UniformEps
  ggsave(filename = "p_UniformEps_AllMethods_MSE_ATT.pdf", p_AllMethods_MSE_ATT_UniformEps,
         width = 26, height = 10,
         units = "cm")

  Figure_UniformEps_BiasSDMSE_ATT <- ggarrange(p_AllMethods_bias_ATT_UniformEps,
                                               p_AllMethods_sd_ATT_UniformEps,
                                               p_AllMethods_MSE_ATT_UniformEps,
                                               nrow = 1, ncol = 3,
                                               widths = c(8, 8, 8, 8), labels = c('(a)', '(b)', '(c)'),
                                               vjust = 1.1, hjust = 0,
                                               common.legend = T, legend = 'right',
                                               font.label = list(size = 23, face = 'plain'))
  Figure_UniformEps_BiasSDMSE_ATT
}
}

# Save .Rdata
save.image(file = "CompareAllMethods_UniformEps_avgEdited.RData")

# Figure Plot: Combine the plot of normal distributed and Uniform distribution----
{
  # Plot of bias+SD+MSE of ATT - Uniform
  Figure_UniformEps_BiasSDMSE_ATT <- ggarrange(p_AllMethods_bias_ATT_UniformEps,
                                               p_AllMethods_sd_ATT_UniformEps,
                                               p_AllMethods_MSE_ATT_UniformEps,
                                               # nrow = 3, ncol = 3,
                                               nrow = 1, ncol = 3,
                                               widths = c(8, 8, 8, 8, 8),
                                               # heights = c(8, 8, 8),
                                               heights = 8,
                                               common.legend = T, legend = 'top',
                                               labels = c('Case II', '', ''),
                                               vjust = 0, hjust = -0.6,
                                               font.label = list(size = 18, face = 'bold'))
  Figure_UniformEps_BiasSDMSE_ATT

  ggsave(filename = "p_UniformEps_BiasSDMSE_ATT.pdf",
         Figure_UniformEps_BiasSDMSE_ATT,
         width = 33, height = 12,
         units = "cm")
}

# -------------------------------------- #
mat_ASMD_CEB_Uniform <- cbind(var_eps = vec_Vareps, mat_ASMD_CEB_Uniform)
colnames(mat_ASMD_CEB_Uniform)[-1] <- paste0(colnames(Ze_Uniform[[1]]), "ASMD")

mat_ASMD_naiveEB_Uniform <- cbind(var_eps = vec_Vareps, mat_ASMD_naiveEB_Uniform)
colnames(mat_ASMD_naiveEB_Uniform)[-1] <- paste0(colnames(Ze_Uniform[[1]]), "ASMD")

mat_ASMD_BCEB_Uniform <- cbind(var_eps = vec_Vareps, mat_ASMD_BCEB_Uniform)
colnames(mat_ASMD_BCEB_Uniform)[-1] <- paste0(colnames(Ze_Uniform[[1]]), "ASMD")

mat_ASMD_HL_Uniform <- cbind(var_eps = vec_Vareps, mat_ASMD_HL_Uniform)
colnames(mat_ASMD_HL_Uniform)[-1] <- paste0(colnames(Ze_Uniform[[1]]), "ASMD")

mat_ASMD_HW_Uniform <- cbind(var_eps = vec_Vareps, mat_ASMD_HW_Uniform)
colnames(mat_ASMD_HW_Uniform)[-1] <- paste0(colnames(Ze_Uniform[[1]]), "ASMD")


round(mat_ASMD_naiveEB_Uniform, 3)
round(mat_ASMD_CEB_Uniform, 3)
round(mat_ASMD_BCEB_Uniform, 3)
round(mat_ASMD_HL_Uniform, 3)
round(mat_ASMD_HW_Uniform, 3)

# Export the data of ASMD into the same Excel with different sheets
library(openxlsx)

list_of_ASMD_UniformME <- list("naiveEB" = round(mat_ASMD_naiveEB_Uniform, 3),
                               "CEB" = round(mat_ASMD_CEB_Uniform, 3),
                               "BCEB" = round(mat_ASMD_BCEB_Uniform, 3),
                               "HL" = round(mat_ASMD_HL_Uniform, 3),
                               "HW" = round(mat_ASMD_HW_Uniform, 3))

openxlsx::write.xlsx(x = list_of_ASMD_UniformME,
                     file = "list_of_ASMD_UniformME.xlsx",
                     row.names = FALSE)

mat_ASMD_naiveEB_Uniform_1 <- cbind(Method = "naiveEB", t(round(mat_ASMD_naiveEB_Uniform, 3)))
mat_ASMD_CEB_Uniform_1 <- cbind(Method = "CEB", t(round(mat_ASMD_CEB_Uniform, 3)))
mat_ASMD_BCEB_Uniform_1 <- cbind(Method = "BCEB", t(round(mat_ASMD_BCEB_Uniform, 3)))
mat_ASMD_HL_Uniform_1 <- cbind(Method = "HL", t(round(mat_ASMD_HL_Uniform, 3)))
mat_ASMD_HW_Uniform_1 <- cbind(Method = "HW", t(round(mat_ASMD_HW_Uniform, 3)))

# 一定要转成dataframe才能成功输出
mat_ASMD_AllMethods <- as.data.frame(rbind(mat_ASMD_naiveEB_Uniform_1,
                                           mat_ASMD_CEB_Uniform_1,
                                           mat_ASMD_BCEB_Uniform_1,
                                           mat_ASMD_HL_Uniform_1,
                                           mat_ASMD_HW_Uniform_1))

openxlsx::write.xlsx(x = mat_ASMD_AllMethods,
                     file = "mat_ASMD_AllMethods_UniformME.xlsx",
                     row.names = TRUE)