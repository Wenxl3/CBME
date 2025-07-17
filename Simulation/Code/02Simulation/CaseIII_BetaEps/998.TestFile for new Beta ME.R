# -*- coding: utf-8 -*-
# ------------------------------------------
# Project     : EBME_R_202306022228
# Title       : 998.TestFile for new Beta ME
# Objective   : Updata Rcode of "Compare all methods with Beta measurement errors"
# Remark      : The measurement error follows a Beta distribution
# Created on  : 2023/06/02 - 22:28
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
  source(file = "./CEB/R/02Simulation/CaseIII_BetaEps/999.ObjectiveGradient_BetaME.R")

  # Variable Information----
{
  ## Sample size
  # N <- 10000
  N <- 2000
  # N <- 5000
  # N <- 500
  # N <- 100
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
  # Beta distributed measurement errors
  measure.error_Beta <- list()
  X_Beta <- list()
  Ze_Beta <- list()
  Ze0_Beta <- list()
  Ze1_Beta <- list()
  Ze1_bar_Beta <- list()
}

  # Simulation Setting----
{
  set.seed(seed = 102)
  ## N_Sim: Number of Monte Carlo experiments
  N_Sim <- 1000
  # N_Sim <- 3000
  # N_Sim <- 500
  # N_Sim <- 20
  # N_Sim <- 5
}

  # Construct data frame for saving results----
{
{
  ### mat_SingleSim: 记录var_eps × rhoZ的某一个组合下的N_Sim次simulation的结果
  #### 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
  #### nrow = N_Sim
  mat_SingleSim_eb_Beta <- matrix(data = 0,
                                  ncol = 7, nrow = N_Sim)
  colnames(mat_SingleSim_eb_Beta) <- c("rhoZ", "Var_eps",
                                       paste0("theta", 1:p_Z),
                                       "ATT")
  mat_SingleSim_ebme_Beta <- mat_SingleSim_eb_Beta
  mat_SingleSim_bc_Beta <- mat_SingleSim_eb_Beta
  mat_SingleSim_hcc_Beta <- mat_SingleSim_eb_Beta
  mat_SingleSim_hyj_Beta <- mat_SingleSim_eb_Beta
  ### mat_avgSim: [每一行]记录var_eps × rhoZ的某一个组合下，500次simulation结果取平均值后的结果
  #### 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
  #### nrow = length(vec_rhoZ) * length(vec_Vareps)
  mat_avgSim_eb_Beta <- matrix(ncol = 7,
                               nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(mat_avgSim_eb_Beta) <- c("rhoZ", "Var_eps",
                                    paste0("theta", 1:p_Z),
                                    "ATT")
  mat_avgSim_ebme_Beta <- mat_avgSim_eb_Beta
  mat_avgSim_bc_Beta <- mat_avgSim_eb_Beta
  mat_avgSim_hcc_Beta <- mat_avgSim_eb_Beta
  mat_avgSim_hyj_Beta <- mat_avgSim_eb_Beta

  ### mat_BiasSdMSE: [每一行]记录var_eps × rhoZ的某一个组合下，500次simulation后估计量的bias, sd, MSE.
  #### 17 cols: rhoZ, Var(eps), bias(theta1, theta2, theta3, theta4, ATT), sd(theta1, theta2, theta3, theta4, ATT), MSE(theta1, theta2, theta3, theta4, ATT)
  #### nrow = length(vec_rhoZ) * length(vec_Vareps)
  mat_BiasSdMSE_eb_Beta <- matrix(ncol = 17,
                                  nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(mat_BiasSdMSE_eb_Beta) <- c("rhoZ", "Var_eps",
                                       paste0("BIAS_", c(paste0("theta", 1:4), "ATT")),
                                       paste0("SD_", c(paste0("theta", 1:4), "ATT")),
                                       paste0("MSE_", c(paste0("theta", 1:4), "ATT")))
  mat_BiasSdMSE_ebme_Beta <- mat_BiasSdMSE_eb_Beta
  mat_BiasSdMSE_bc_Beta <- mat_BiasSdMSE_eb_Beta
  mat_BiasSdMSE_hcc_Beta <- mat_BiasSdMSE_eb_Beta
  mat_BiasSdMSE_hyj_Beta <- mat_BiasSdMSE_eb_Beta

  ## Iteration numbers and indicator pool
  ### EBME-BFGS
  iterbfgsind.r_Beta <- rep(TRUE, N_Sim)
  ### BC - approximate bias corrector
  iterbcind.r_Beta <- iterbfgsind.r_Beta
  ### HCC
  iterhccind.r_Beta <- iterbfgsind.r_Beta
  ### HYJ
  iterhyjind.r_Beta <- iterbfgsind.r_Beta

  ## Rate of model convergence: Calculate the rate of convergence in each combination of rho_Z and var_eps
  ### 3 cols: rhoZ, Var(eps), average value of the rate of convergence
  ### nrow <- length(vec_rhoZ) * length(vec_Vareps)
  Rc_ebme_Beta <- matrix(data = 0,
                         ncol = 3,
                         nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(Rc_ebme_Beta) <- c("rhoZ", "Var_eps", "Rc")
  Rc_bc_Beta <- Rc_ebme_Beta
  Rc_hcc_Beta <- Rc_ebme_Beta
  Rc_hyj_Beta <- Rc_ebme_Beta

  # matrix markdown the ASMD
  mat_ASMD_CEB_Beta <- NULL
  mat_ASMD_naiveEB_Beta <- NULL
  mat_ASMD_BCEB_Beta <- NULL
  mat_ASMD_HL_Beta <- NULL
  mat_ASMD_HW_Beta <- NULL
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
    # var_eps <- 0.5

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
      ## Beta distributed measurement errors
      mat_SingleSim_eb_Beta[i, "rhoZ"] <- rhoZ
      mat_SingleSim_ebme_Beta[i, "rhoZ"] <- rhoZ
      mat_SingleSim_bc_Beta[i, "rhoZ"] <- rhoZ
      mat_SingleSim_hcc_Beta[i, "rhoZ"] <- rhoZ
      mat_SingleSim_hyj_Beta[i, "rhoZ"] <- rhoZ
    }

      # Mark the Var_eps
    {
      ## Beta distributed measurement errors
      mat_SingleSim_eb_Beta[i, "Var_eps"] <- var_eps
      mat_SingleSim_ebme_Beta[i, "Var_eps"] <- var_eps
      mat_SingleSim_bc_Beta[i, "Var_eps"] <- var_eps
      mat_SingleSim_hcc_Beta[i, "Var_eps"] <- var_eps
      mat_SingleSim_hyj_Beta[i, "Var_eps"] <- var_eps
    }

      # True covariates
    {
      Z <- mvrnorm(n = N, mu = c(Mean_X, Mean_U), Sigma = Sigma_Z)
      trueX <- Z[, 1:p_X]
      U <- Z[, (p_X + 1):p_Z]
    }

      # Construct replicated covariate dataset
    {
      # Beta(3,1) measurement error
      shape_1 <- 3
      shape_2 <- 1
      # mean value of the uniform distribution
      mean_val_Beta <- shape_1 / (shape_1 + shape_2)
      # variance value of the uniform distribution
      var_val_Beta <- shape_1 * shape_2 / ((shape_1 + shape_2)^2 * (shape_1 + shape_2 + 1))

      for (k in 1:r) {
        if (var_eps == 0) {
          measure.error_Beta[[k]] <- rep(0, p_X * N)
        }else {
          measure.error_Beta[[k]] <- rbeta(n = p_X * N, shape1 = shape_1, shape2 = shape_2)

          # Adjust the expectation of measurement error to 0
          measure.error_Beta[[k]] <- measure.error_Beta[[k]] - mean_val_Beta
          # Adjust the variance of measurement error to prespecified variance
          measure.error_Beta[[k]] <- measure.error_Beta[[k]] * sqrt(var_eps) / sqrt(var_val_Beta)
        }

        measure.error_Beta[[k]] <- matrix(data = measure.error_Beta[[k]], nrow = N, ncol = p_X)

        X_Beta[[k]] <- trueX + measure.error_Beta[[k]]
        Ze_Beta[[k]] <- data.frame(cbind(X_Beta[[k]], U))
        colnames(Ze_Beta[[k]]) <- c("X1", "X2", "U3", "U4")
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
      Ze1_bar_al_Beta <- rep(0, p_Z)
      for (k in 1:r) {
        Ze0_Beta[[k]] <- as.matrix(Ze_Beta[[k]][Index_col,])
        Ze1_Beta[[k]] <- as.matrix(Ze_Beta[[k]][Index_trt,])
        Ze1_bar_Beta[[k]] <- colMeans(Ze1_Beta[[k]])
        Ze1_bar_al_Beta <- Ze1_bar_al_Beta + Ze1_bar_Beta[[k]]
        Ze_Beta[[k]] <- as.matrix(Ze_Beta[[k]])
      }
      Ze1_bar_al_Beta <- Ze1_bar_al_Beta / r
    }

      # ZQY-EB: native Entropy Balancing
    {
      opt.out.zqy.r_Beta <- optim(par = rep(0, (ncol(Z))),
                                  fn = objective.zqy.r_Beta,
                                  gr = gradient.zqy.r_Beta,
                                  method = "BFGS",
                                  control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200))
      theta0_Beta <- opt.out.zqy.r_Beta$par
      # round(theta0_Beta, 3)
    }


      # EBME - BFGS
    {
      # theta_CEB - solve the objective funtion
      output.CEB.solveObj_Beta <- tryCatch(optim(par = rep(0, p_Z),
                                                 fn = objective.ebme.r_Beta,
                                                 gr = gradient.ebme.r_Beta,
                                                 method = "BFGS",
                                                 control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200)),
                                           error = function(e) {
                                             print("Fail to solve objective of CEB, We directly use naive EB's results for subsequent comparison.")
                                             # 此时无法求解CEB，直接将solve CEB Objective的结果赋为FALSE
                                             return(FALSE)
                                           })
      # 如果optim无法求解，output.CEB.solveObj_Beta会是一个logical变量；否则，如果optim能求解，output.CEB.solveObj_Beta是一个list变量
      if (class(output.CEB.solveObj_Beta) != "list" |
        is.na(output.CEB.solveObj_Beta$value) |
        is.nan(output.CEB.solveObj_Beta$value)) {
        # 此时无法求解CEB，将naive EB的结果赋给CEB
        output.CEB.solveObj_Beta <- opt.out.zqy.r_Beta
      }
      theta_CEB_solveObj <- output.CEB.solveObj_Beta$par
      ObjValue_CEB_solveObj <- output.CEB.solveObj_Beta$value

      # theta_CEB - solve the gradient funtion
      output.CEB.solveGradient_Beta <- tryCatch(multiroot(f = gradient.ebme.r_Beta,
                                                          start = rep(0, p_Z),
                                                          rtol = .Machine$double.eps, maxiter = 200),
                                                error = function(e) {
                                                  print("Fail to solve gradient of CEB, We directly use naive EB's results for subsequent comparison.")
                                                  # 此时无法求解CEB，直接将solve CEB Objective的结果赋为FALSE
                                                  return(FALSE)
                                                })
      # 如果multiroot无法求解，output.CEB.solveGradient会是一个logical变量；否则，如果multiroot能求解，output.CEB.solveGradient是一个list变量
      if (class(output.CEB.solveGradient_Beta) != "list") {
        # 此时求解CEB gradient 不收敛，直接赋值9999
        ObjValue_CEB_solveGradient <- 9999
      }else {
        theta_CEB_solveGradient <- output.CEB.solveGradient_Beta$root
        ObjValue_CEB_solveGradient <- objective.ebme.r_Beta(theta = theta_CEB_solveGradient)
        if (is.nan(ObjValue_CEB_solveGradient) | is.na(ObjValue_CEB_solveGradient)) { ObjValue_CEB_solveGradient <- 9999 }
      }

      if (ObjValue_CEB_solveObj <= ObjValue_CEB_solveGradient) {
        theta_CEB_Beta <- theta_CEB_solveObj
        iterbfgsind.r_Beta[i] <- (max(abs(gradient.ebme.r_Beta(theta_CEB_solveObj))) < 1e-04)
        ## double check
        if (is.na(iterbfgsind.r_Beta[i]) | is.nan(iterbfgsind.r_Beta[i])) {
          iterbfgsind.r_Beta[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }else {
        theta_CEB_Beta <- theta_CEB_solveGradient
        iterbfgsind.r_Beta[i] <- (max(abs(output.CEB.solveGradient_Beta$f.root)) < sqrt(.Machine$double.eps)) & (output.CEB.solveGradient_Beta$estim.precis != "NaN")
        ## double check
        if (is.na(iterbfgsind.r_Beta[i]) | is.nan(iterbfgsind.r_Beta[i])) {
          iterbfgsind.r_Beta[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }
    }

      # BCEB
    {
      theta_bc_Beta <- drop(solve(hessian.zqy.r_Beta(theta0_Beta) - Sigma) %*%
                              hessian.zqy.r_Beta(theta0_Beta) %*%
                              theta0_Beta)
      ### 判断BC是否无法求解：（默认【收敛】）存在以下情况任意其一都被判为【无解/不收敛】
      iterbcind.r_Beta[i] <- !(NaN %in% theta_bc_Beta | (sum(is.infinite(theta_bc_Beta)) > 0))
    }

      # HCC(2012)
    {
      opt.out.hcc.r_Beta <- multiroot(f = gradient.hcc.r_Beta, start = rep(0, p_Z),
                                      rtol = .Machine$double.eps, maxiter = 200)
      iterhccind.r_Beta[i] <- (max(gradient.hcc.r_Beta(opt.out.hcc.r_Beta$root)) < sqrt(.Machine$double.eps)) & (opt.out.hcc.r_Beta$estim.precis != "NaN")
    }

      # HYJ(2012)
    {
      opt.out.hyj.r_Beta <- multiroot(f = gradient.hyj.r_Beta, start = rep(0, p_Z),
                                      rtol = .Machine$double.eps, maxiter = 200)
      iterhyjind.r_Beta[i] <- (max(gradient.hyj.r_Beta(opt.out.hyj.r_Beta$root)) < sqrt(.Machine$double.eps)) & (opt.out.hyj.r_Beta$estim.precis != "NaN")
      ### double check
      if (is.na(iterhyjind.r_Beta[i])) {
        iterhyjind.r_Beta[i] <- FALSE
        print("NA clear!")
      }
    }

      # ATT estimation with \hat{\theta}
    {
      # ZQY-EB: naive Entropy Balancing
    {
      # par_zqy_Beta <- as.matrix(opt.out.zqy.r_Beta$par, ncol = 1)
      # weight_zqy_Beta <- exp(Ze0_Beta[[1]] %*% par_zqy_Beta)
      weight_zqy_Beta <- exp(Ze0_Beta[[1]] %*% theta0_Beta)
      weight_zqy_Beta <- weight_zqy_Beta / sum(weight_zqy_Beta)
      ATT_zqy_Beta <- mean(Y1) - t(weight_zqy_Beta) %*% Y0
      mat_SingleSim_eb_Beta[i, "ATT"] <- mat_SingleSim_eb_Beta[i, "ATT"] + ATT_zqy_Beta
    }

      # CEB / EBME - bfgs C2: the proposed method
    {
      # par_bfgs_Beta <- as.matrix(theta_CEB_Beta, ncol = 1)
      # weight_bfgs_Beta <- exp(Ze0_Beta[[1]] %*% par_bfgs_Beta)
      weight_bfgs_Beta <- exp(Ze0_Beta[[1]] %*% theta_CEB_Beta)
      weight_bfgs_Beta <- weight_bfgs_Beta / sum(weight_bfgs_Beta)
      ATT_bfgs_Beta <- mean(Y1) - t(weight_bfgs_Beta) %*% Y0
      mat_SingleSim_ebme_Beta[i, "ATT"] <- mat_SingleSim_ebme_Beta[i, "ATT"] + ATT_bfgs_Beta
    }

      # BCEB
    {
      # par_bc_Beta <- as.matrix(theta_bc_Beta, ncol = 1)
      # weight_bc_Beta <- exp(Ze0_Beta[[1]] %*% par_bc_Beta)
      weight_bc_Beta <- exp(Ze0_Beta[[1]] %*% theta_bc_Beta)
      weight_bc_Beta <- weight_bc_Beta / sum(weight_bc_Beta)
      ATT_bc_Beta <- mean(Y1) - t(weight_bc_Beta) %*% Y0
      mat_SingleSim_bc_Beta[i, "ATT"] <- mat_SingleSim_bc_Beta[i, "ATT"] + ATT_bc_Beta
    }

      # HCC
    {
      par_hcc_Beta <- as.matrix(opt.out.hcc.r_Beta$root, ncol = 1)
      weight_hcc_Beta <- matrix(rep(0, nrow(weight_bc_Beta)))
      for (k in 1:r) {
        weight_hcc_Beta <- weight_hcc_Beta + exp(Ze0_Beta[[k]] %*% par_hcc_Beta)
      }
      weight_hcc_Beta <- weight_hcc_Beta / sum(weight_hcc_Beta)
      ATT_hcc_Beta <- mean(Y1) - t(weight_hcc_Beta) %*% Y0
      mat_SingleSim_hcc_Beta[i, "ATT"] <- mat_SingleSim_hcc_Beta[i, "ATT"] + ATT_hcc_Beta
    }

      # HYJ
    {
      par_hyj_Beta <- as.matrix(opt.out.hyj.r_Beta$root, ncol = 1)
      weight_hyj_Beta <- matrix(rep(0, nrow(weight_bc_Beta)))
      for (k in 1:r) {
        weight_hyj_Beta <- weight_hyj_Beta + exp(Ze0_Beta[[k]] %*% par_hyj_Beta)
      }
      weight_hyj_Beta <- weight_hyj_Beta / sum(weight_hyj_Beta)
      ATT_hyj_Beta <- mean(Y1) - t(weight_hyj_Beta) %*% Y0
      ### double check
      if (is.na(ATT_hyj_Beta)) {
        iterhyjind.r_Beta[i] <- FALSE
        print("ATT of HW - NA clear!")
      }
      mat_SingleSim_hyj_Beta[i, "ATT"] <- mat_SingleSim_hyj_Beta[i, "ATT"] + ATT_hyj_Beta
    }

    }

      # theta coefficient
    {
      ## Beta distributed measurement errors
      mat_SingleSim_eb_Beta[i, paste0("theta", 1:p_Z)] <- theta0_Beta
      mat_SingleSim_ebme_Beta[i, paste0("theta", 1:p_Z)] <- theta_CEB_Beta
      mat_SingleSim_bc_Beta[i, paste0("theta", 1:p_Z)] <- theta_bc_Beta
      mat_SingleSim_hcc_Beta[i, paste0("theta", 1:p_Z)] <- opt.out.hcc.r_Beta$root
      mat_SingleSim_hyj_Beta[i, paste0("theta", 1:p_Z)] <- opt.out.hyj.r_Beta$root
    }

      # ASMD calculation ----
    {
      # mat_SingleSim_ASMD_CEB
      mat_SingleSim_ASMD_CEB[i,] <- unlist(ASMDCal(weight_bfgs_Beta))
      mat_SingleSim_ASMD_naiveEB[i,] <- unlist(ASMDCal(weight_zqy_Beta))
      mat_SingleSim_ASMD_BCEB[i,] <- unlist(ASMDCal(weight_bc_Beta))
      mat_SingleSim_ASMD_HL[i,] <- unlist(ASMDCal(weight_hcc_Beta))
      mat_SingleSim_ASMD_HW[i,] <- unlist(ASMDCal(weight_hyj_Beta))
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
    Rc_ebme_Beta[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterbfgsind.r_Beta) / N_Sim)
    ## BC
    Rc_bc_Beta[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterbcind.r_Beta) / N_Sim)
    ## HCC
    Rc_hcc_Beta[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterhccind.r_Beta) / N_Sim)
    ## HYJ
    Rc_hyj_Beta[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterhyjind.r_Beta) / N_Sim)
  }

    # mat_avgSim赋值 - Convergent results
  {
    mat_avgSim_eb_Beta[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_eb_Beta)
    mat_avgSim_ebme_Beta[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_ebme_Beta[iterbfgsind.r_Beta,])
    mat_avgSim_bc_Beta[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_bc_Beta[iterbcind.r_Beta,])
    mat_avgSim_hcc_Beta[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_hcc_Beta[iterhccind.r_Beta,])
    mat_avgSim_hyj_Beta[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_hyj_Beta[iterhyjind.r_Beta,])
  }

    # mat_BiasSDMSE赋值 - Convergent results
    trueVal <- c(True.theta[-1], ATE)
  {
    ## rhoZ赋值
    mat_BiasSdMSE_eb_Beta[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_ebme_Beta[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_bc_Beta[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_hcc_Beta[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_hyj_Beta[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    ## Var_eps
    mat_BiasSdMSE_eb_Beta[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_ebme_Beta[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_bc_Beta[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_hcc_Beta[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_hyj_Beta[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    ## bias
    ColnamesinNeed_avgmat <- c(paste0("theta", 1:p_Z), "ATT")
    ColnamesinNeed_Bias <- paste0("BIAS_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Beta[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_eb_Beta[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_ebme_Beta[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_ebme_Beta[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_bc_Beta[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_bc_Beta[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_hcc_Beta[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_hcc_Beta[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_hyj_Beta[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_hyj_Beta[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    ## sd
    ColnamesinNeed_SD <- paste0("SD_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Beta[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_eb_Beta[, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_ebme_Beta[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_ebme_Beta[iterbfgsind.r_Beta, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd, na.rm = TRUE)
    mat_BiasSdMSE_bc_Beta[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_bc_Beta[iterbcind.r_Beta, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd, na.rm = TRUE)
    mat_BiasSdMSE_hcc_Beta[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_hcc_Beta[iterhccind.r_Beta, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_hyj_Beta[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_hyj_Beta[iterhyjind.r_Beta, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    ## MSE
    ColnamesinNeed_MSE <- paste0("MSE_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Beta[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_eb_Beta[, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_ebme_Beta[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_ebme_Beta[iterbfgsind.r_Beta, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_bc_Beta[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_bc_Beta[iterbcind.r_Beta, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_hcc_Beta[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_hcc_Beta[iterhccind.r_Beta, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_hyj_Beta[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_hyj_Beta[iterhyjind.r_Beta, ColnamesinNeed_avgmat])
  }

    # ASMD
  {
    # mat_ASMD_CEB_Beta
    mat_ASMD_naiveEB_Beta <- rbind(mat_ASMD_naiveEB_Beta,
                                   colMeans(mat_SingleSim_ASMD_naiveEB, na.rm = TRUE))
    mat_ASMD_CEB_Beta <- rbind(mat_ASMD_CEB_Beta,
                               colMeans(mat_SingleSim_ASMD_CEB, na.rm = TRUE))
    mat_ASMD_BCEB_Beta <- rbind(mat_ASMD_BCEB_Beta,
                                colMeans(mat_SingleSim_ASMD_BCEB, na.rm = TRUE))
    mat_ASMD_HL_Beta <- rbind(mat_ASMD_HL_Beta,
                              colMeans(mat_SingleSim_ASMD_HL, na.rm = TRUE))
    mat_ASMD_HW_Beta <- rbind(mat_ASMD_HW_Beta,
                              colMeans(mat_SingleSim_ASMD_HW, na.rm = TRUE))
    # round(mat_ASMD_CEB_Beta, 3)
  }


    # 清空记录上一个var_eps × rhoZ组合的500次simulation结果的数据集
    Count_Row <- 0
  {
  {
    # mat_SingleSim: 记录var_eps × rhoZ的某一个组合下的N_Sim次simulation的结果
    ## 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
    ## nrow = N_Sim
    mat_SingleSim_eb_Beta <- matrix(data = 0,
                                    ncol = 7, nrow = N_Sim)
    colnames(mat_SingleSim_eb_Beta) <- c("rhoZ", "Var_eps",
                                         paste0("theta", 1:p_Z), "ATT")
    mat_SingleSim_ebme_Beta <- mat_SingleSim_eb_Beta
    mat_SingleSim_bc_Beta <- mat_SingleSim_eb_Beta
    mat_SingleSim_hcc_Beta <- mat_SingleSim_eb_Beta
    mat_SingleSim_hyj_Beta <- mat_SingleSim_eb_Beta

    # 清空记录上一个var_eps × rhoZ组合下500次simulation的收敛结果标记
    ## EBME-BFGS
    iterbfgsind.r_Beta <- rep(TRUE, N_Sim)
    ## BC - approximate bias corrector
    iterbcind.r_Beta <- iterbfgsind.r_Beta
    ## HCC
    iterhccind.r_Beta <- iterbfgsind.r_Beta
    ## HYJ
    iterhyjind.r_Beta <- iterbfgsind.r_Beta


    # 矩阵清空为下一轮simulation做准备
    mat_SingleSim_ASMD_CEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_naiveEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_BCEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HL <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HW <- matrix(0, nrow = N_Sim, ncol = p_Z)
  }
  }
  }
}


# Save .Rdata
save.image(file = "CompareAllMethods_BetaEps_Original.RData")

{
  # Table output----
{
  Rc_ebme_Beta <- cbind(Rc_ebme_Beta, Method = "CEB")
  Rc_bc_Beta <- cbind(Rc_bc_Beta, Method = "BCEB")
  Rc_hcc_Beta <- cbind(Rc_hcc_Beta, Method = "HL")
  Rc_hyj_Beta <- cbind(Rc_hyj_Beta, Method = "HW")
  mat_Rc_Beta <- rbind(Rc_ebme_Beta,
                       Rc_bc_Beta,
                       Rc_hcc_Beta,
                       Rc_hyj_Beta)
  csvName_mat_Rc_Beta <- paste0("mat_BetaEps_ConvergenceRate_N", N,
                                "_Sim", N_Sim, ".csv")
  write.table(mat_Rc_Beta, file = csvName_mat_Rc_Beta,
              sep = ",", row.names = F)

  ## mat_avgSim
  mat_avgSim_eb_Beta <- cbind(mat_avgSim_eb_Beta, Method = "Naive")
  mat_avgSim_ebme_Beta <- cbind(mat_avgSim_ebme_Beta, Method = "CEB")
  mat_avgSim_bc_Beta <- cbind(mat_avgSim_bc_Beta, Method = "BCEB")
  mat_avgSim_hcc_Beta <- cbind(mat_avgSim_hcc_Beta, Method = "HL")
  mat_avgSim_hyj_Beta <- cbind(mat_avgSim_hyj_Beta, Method = "HW")
  mat_avgSim_Beta <- rbind(mat_avgSim_eb_Beta,
                           mat_avgSim_ebme_Beta,
                           mat_avgSim_bc_Beta,
                           mat_avgSim_hcc_Beta,
                           mat_avgSim_hyj_Beta)
  csvName_mat_avgSim <- paste0("mat_BetaEps_avgSim_theta_ATT_N", N, "_Sim", N_Sim,
                               ".csv")
  write.table(mat_avgSim_Beta, file = csvName_mat_avgSim,
              sep = ",", row.names = F)

  ## mat_BiasSdMSE
  mat_BiasSdMSE_eb_Beta <- cbind(mat_BiasSdMSE_eb_Beta, Method = "Naive")
  mat_BiasSdMSE_ebme_Beta <- cbind(mat_BiasSdMSE_ebme_Beta, Method = "CEB")
  mat_BiasSdMSE_bc_Beta <- cbind(mat_BiasSdMSE_bc_Beta, Method = "BCEB")
  mat_BiasSdMSE_hcc_Beta <- cbind(mat_BiasSdMSE_hcc_Beta, Method = "HL")
  mat_BiasSdMSE_hyj_Beta <- cbind(mat_BiasSdMSE_hyj_Beta, Method = "HW")
  mat_BiasSdMSE_Beta <- rbind(mat_BiasSdMSE_eb_Beta,
                              mat_BiasSdMSE_ebme_Beta,
                              mat_BiasSdMSE_bc_Beta,
                              mat_BiasSdMSE_hcc_Beta,
                              mat_BiasSdMSE_hyj_Beta)
  csvName_mat_BiasSdMSE <- paste0("mat_BetaEps_BiasSdMSE_theta_ATT_N", N, "_Sim", N_Sim,
                                  ".csv")
  write.table(mat_BiasSdMSE_Beta, file = csvName_mat_BiasSdMSE,
              sep = ",", row.names = F)
}

  # Plots----
{
  # Plot of convergenc rate of EBME, BC, HCC, HYJ
  Rcdf_Beta <- as.data.frame(mat_Rc_Beta)
  Rcdf_Beta$Method <- factor(Rcdf_Beta$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  p_Rc_BetaEps <- ggplot(Rcdf_Beta) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = Rc, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = Rc, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Rate of Convergence", title = "Measurement Error with a Beta Distribution") +
    theme_bw()
  p_Rc_BetaEps
  ggsave(filename = "p_BetaEps_Rc.pdf", p_Rc_BetaEps,
         width = 14, height = 13,
         units = "cm")

  ## Seperate bias, sd, MSE dataframe
  mat_bias_Beta <- mat_BiasSdMSE_Beta[, c("rhoZ", "Var_eps", paste0("BIAS_", c(paste0("theta", 1:4), "ATT")), "Method")]
  mat_sd_Beta <- mat_BiasSdMSE_Beta[, c("rhoZ", "Var_eps", paste0("SD_", c(paste0("theta", 1:4), "ATT")), "Method")]
  mat_MSE_Beta <- mat_BiasSdMSE_Beta[, c("rhoZ", "Var_eps", paste0("MSE_", c(paste0("theta", 1:4), "ATT")), "Method")]
  ## bias
  biasdf_Beta <- melt(as.data.frame(mat_bias_Beta),
                      id.vars = c("rhoZ", "Var_eps", "Method"))
  biasdf_Beta$value <- as.numeric(biasdf_Beta$value)
  which(is.na(biasdf_Beta$value))
  # biasdf_Beta$Method <- factor(biasdf_Beta$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  biasdf_Beta$Method <- factor(biasdf_Beta$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### bias_theta
  biasdf_theta_Beta <- biasdf_Beta[biasdf_Beta$variable != "BIAS_ATT",]
  p_AllMethods_bias_theta_BetaEps <- ggplot(biasdf_theta_Beta) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) + # 行按variable进行分区，列按rhoZ分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Bias", title = "True Value: theta1 = -1, theta2 = 0.5, theta2 = -0.25, theta4 = -0.1") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    # ylim(-0.8, 0.6) +
    # scale_y_continuous(breaks = seq(-0.8, 0.6, by = 0.1)) +
    theme_bw()
  p_AllMethods_bias_theta_BetaEps
  ggsave(filename = "p_BetaEps_AllMethods_bias_theta.pdf", p_AllMethods_bias_theta_BetaEps,
         width = 26, height = 26,
         units = "cm")
  ### bias_ATT
  biasdf_ATT_Beta <- biasdf_Beta[biasdf_Beta$variable == "BIAS_ATT",]
  p_AllMethods_bias_ATT_BetaEps <- ggplot(biasdf_ATT_Beta) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    # labs(x = "Variance of Measurement Error", y = "Bias", title = "Bias of ATT with Beta Measurement Errors") +
    labs(x = "Variance of Measurement Error", y = "Bias") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    # ylim(-9, 2) +
    # scale_y_continuous(breaks = seq(-9, 2, by = 0.5)) +
    theme_bw()
  p_AllMethods_bias_ATT_BetaEps
  ggsave(filename = "p_BetaEps_AllMethods_bias_ATT.pdf", p_AllMethods_bias_ATT_BetaEps,
         width = 26, height = 10,
         units = "cm")
  ## sd
  sddf_Beta <- melt(as.data.frame(mat_sd_Beta),
                    id.vars = c("rhoZ", "Var_eps", "Method"))
  sddf_Beta$value <- as.numeric(sddf_Beta$value)
  which(is.na(sddf_Beta$value))
  # sddf_Beta$Method <- factor(sddf_Beta$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  sddf_Beta$Method <- factor(sddf_Beta$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### sd_theta
  sddf_theta_Beta <- sddf_Beta[sddf_Beta$variable != "SD_ATT",]
  p_AllMethods_sd_theta_BetaEps <- ggplot(sddf_theta_Beta) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) + # 行按variable进行分区，列按rhoZ分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Standard Deviation", title = "SD of Coefficients with Beta Distributed Measurement Error") +
    theme_bw()
  p_AllMethods_sd_theta_BetaEps
  ggsave(filename = "p_BetaEps_AllMethods_sd_theta.pdf", p_AllMethods_sd_theta_BetaEps,
         width = 26, height = 26,
         units = "cm")
  ### sd_ATT
  sddf_ATT_Beta <- sddf_Beta[sddf_Beta$variable == "SD_ATT",]
  p_AllMethods_sd_ATT_BetaEps <- ggplot(sddf_ATT_Beta) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    # scale_y_continuous(limits = c(0, 2.25), breaks = seq(from = 0, to = 2.25, by = 0.5)) +
    scale_y_continuous(limits = c(0, max(sddf_ATT_Beta$value)), breaks = seq(from = 0, to = ceiling(max(sddf_ATT_Beta$value)), by = 0.5)) +
    # labs(x = "Variance of Measurement Error", y = "Standard Deviation", title = "SD of ATT with Beta Measurement Errors") +
    labs(x = "Variance of Measurement Error", y = "Standard Deviation") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_sd_ATT_BetaEps
  ggsave(filename = "p_BetaEps_AllMethods_sd_ATT.pdf", p_AllMethods_sd_ATT_BetaEps,
         width = 26, height = 10,
         units = "cm")
  ## MSE
  MSEdf_Beta <- melt(as.data.frame(mat_MSE_Beta),
                     id.vars = c("rhoZ", "Var_eps", "Method"))
  MSEdf_Beta$value <- as.numeric(MSEdf_Beta$value)
  which(is.na(MSEdf_Beta$value))
  # MSEdf_Beta$Method <- factor(MSEdf_Beta$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  MSEdf_Beta$Method <- factor(MSEdf_Beta$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### MSE_theta
  MSEdf_theta_Beta <- MSEdf_Beta[MSEdf_Beta$variable != "MSE_ATT",]
  p_AllMethods_MSE_theta_BetaEps <- ggplot(MSEdf_theta_Beta) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Mean Squared Error", title = "MSE of Coefficients with Beta Distributed Measurement Error") +
    theme_bw()
  p_AllMethods_MSE_theta_BetaEps
  ggsave(filename = "p_BetaEps_AllMethods_MSE_theta.pdf", p_AllMethods_MSE_theta_BetaEps,
         width = 26, height = 26,
         units = "cm")
  ### MSE_ATT
  MSEdf_ATT_Beta <- MSEdf_Beta[MSEdf_Beta$variable == "MSE_ATT",]
  p_AllMethods_MSE_ATT_BetaEps <- ggplot(MSEdf_ATT_Beta) +
    # facet_wrap(~rhoZ) +  # 行按rhoZ进行分区
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    # labs(x = "Variance of Measurement Error", y = "Mean Squared Error", title = "MSE of ATT with Beta Measurement Errors") +
    labs(x = "Variance of Measurement Error", y = "Mean Squared Error") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_MSE_ATT_BetaEps
  ggsave(filename = "p_BetaEps_AllMethods_MSE_ATT.pdf", p_AllMethods_MSE_ATT_BetaEps,
         width = 26, height = 10,
         units = "cm")

  Figure_BetaEps_BiasSDMSE_ATT <- ggarrange(p_AllMethods_bias_ATT_BetaEps,
                                            p_AllMethods_sd_ATT_BetaEps,
                                            p_AllMethods_MSE_ATT_BetaEps,
                                            nrow = 1, ncol = 3,
                                            widths = c(8, 8, 8, 8), labels = c('(a)', '(b)', '(c)'),
                                            vjust = 1.1, hjust = 0,
                                            common.legend = T, legend = 'right',
                                            font.label = list(size = 23, face = 'plain'))
  Figure_BetaEps_BiasSDMSE_ATT
}
}

# Save .Rdata
save.image(file = "CompareAllMethods_BetaEps_avgEdited.RData")

# Figure Plot: Combine the plot of Beta distributed and Beta distribution----
{
  # Plot of bias+SD+MSE of ATT - Beta
  Figure_BetaEps_BiasSDMSE_ATT <- ggarrange(p_AllMethods_bias_ATT_BetaEps,
                                            p_AllMethods_sd_ATT_BetaEps,
                                            p_AllMethods_MSE_ATT_BetaEps,
                                            # nrow = 3, ncol = 3,
                                            nrow = 1, ncol = 3,
                                            widths = c(8, 8, 8, 8, 8),
                                            # heights = c(8, 8, 8),
                                            heights = 8,
                                            common.legend = T, legend = 'top',
                                            labels = c('Case III', '', ''),
                                            vjust = 0, hjust = -0.6,
                                            font.label = list(size = 18, face = 'bold'))
  Figure_BetaEps_BiasSDMSE_ATT

  ggsave(filename = "p_BetaEps_BiasSDMSE_ATT.pdf",
         Figure_BetaEps_BiasSDMSE_ATT,
         width = 33, height = 12,
         units = "cm")
}


# -------------------------------------- #
mat_ASMD_CEB_Beta <- cbind(var_eps = vec_Vareps, mat_ASMD_CEB_Beta)
colnames(mat_ASMD_CEB_Beta)[-1] <- paste0(colnames(Ze_Beta[[1]]), "ASMD")

mat_ASMD_naiveEB_Beta <- cbind(var_eps = vec_Vareps, mat_ASMD_naiveEB_Beta)
colnames(mat_ASMD_naiveEB_Beta)[-1] <- paste0(colnames(Ze_Beta[[1]]), "ASMD")

mat_ASMD_BCEB_Beta <- cbind(var_eps = vec_Vareps, mat_ASMD_BCEB_Beta)
colnames(mat_ASMD_BCEB_Beta)[-1] <- paste0(colnames(Ze_Beta[[1]]), "ASMD")

mat_ASMD_HL_Beta <- cbind(var_eps = vec_Vareps, mat_ASMD_HL_Beta)
colnames(mat_ASMD_HL_Beta)[-1] <- paste0(colnames(Ze_Beta[[1]]), "ASMD")

mat_ASMD_HW_Beta <- cbind(var_eps = vec_Vareps, mat_ASMD_HW_Beta)
colnames(mat_ASMD_HW_Beta)[-1] <- paste0(colnames(Ze_Beta[[1]]), "ASMD")


round(mat_ASMD_naiveEB_Beta, 3)
round(mat_ASMD_CEB_Beta, 3)
round(mat_ASMD_BCEB_Beta, 3)
round(mat_ASMD_HL_Beta, 3)
round(mat_ASMD_HW_Beta, 3)

# Export the data of ASMD into the same Excel with different sheets
library(openxlsx)

list_of_ASMD_BetaME <- list("naiveEB" = round(mat_ASMD_naiveEB_Beta, 3),
                            "CEB" = round(mat_ASMD_CEB_Beta, 3),
                            "BCEB" = round(mat_ASMD_BCEB_Beta, 3),
                            "HL" = round(mat_ASMD_HL_Beta, 3),
                            "HW" = round(mat_ASMD_HW_Beta, 3))

openxlsx::write.xlsx(x = list_of_ASMD_BetaME,
                     file = "list_of_ASMD_BetaME.xlsx",
                     row.names = FALSE)

mat_ASMD_naiveEB_Beta_1 <- cbind(Method = "naiveEB", t(round(mat_ASMD_naiveEB_Beta, 3)))
mat_ASMD_CEB_Beta_1 <- cbind(Method = "CEB", t(round(mat_ASMD_CEB_Beta, 3)))
mat_ASMD_BCEB_Beta_1 <- cbind(Method = "BCEB", t(round(mat_ASMD_BCEB_Beta, 3)))
mat_ASMD_HL_Beta_1 <- cbind(Method = "HL", t(round(mat_ASMD_HL_Beta, 3)))
mat_ASMD_HW_Beta_1 <- cbind(Method = "HW", t(round(mat_ASMD_HW_Beta, 3)))

# 一定要转成dataframe才能成功输出
mat_ASMD_AllMethods <- as.data.frame(rbind(mat_ASMD_naiveEB_Beta_1,
                                           mat_ASMD_CEB_Beta_1,
                                           mat_ASMD_BCEB_Beta_1,
                                           mat_ASMD_HL_Beta_1,
                                           mat_ASMD_HW_Beta_1))

openxlsx::write.xlsx(x = mat_ASMD_AllMethods,
                     file = "mat_ASMD_AllMethods_BetaME.xlsx",
                     row.names = TRUE)