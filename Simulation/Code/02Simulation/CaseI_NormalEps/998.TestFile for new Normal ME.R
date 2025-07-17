# -*- coding: utf-8 -*-
# ------------------------------------------
# Title       : 998.TestFile for new Normal ME
# Remark      : The measurement error follows N(0,Sigma).
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
  source(file = "./CEB/R/02Simulation/CaseI_NormalEps/999.ObjectiveGradient_NormalME.R")

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
  ## Normal distributed measurement errors
  measure.error_Normal <- list()
  X_Normal <- list()
  Ze_Normal <- list()
  Ze0_Normal <- list()
  Ze1_Normal <- list()
  Ze1_bar_Normal <- list()
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
  ### mat_SingleSim: 记录var_eps × rhoZ的某一个组合下的N_Sim次simulation的结果
  #### 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
  #### nrow = N_Sim
  mat_SingleSim_eb_Normal <- matrix(data = 0,
                                    ncol = 7, nrow = N_Sim)
  colnames(mat_SingleSim_eb_Normal) <- c("rhoZ", "Var_eps",
                                         paste0("theta", 1:p_Z),
                                         "ATT")
  mat_SingleSim_ebme_Normal <- mat_SingleSim_eb_Normal
  mat_SingleSim_bc_Normal <- mat_SingleSim_eb_Normal
  mat_SingleSim_hcc_Normal <- mat_SingleSim_eb_Normal
  mat_SingleSim_hyj_Normal <- mat_SingleSim_eb_Normal

  ### mat_avgSim:
  #### 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
  #### nrow = length(vec_rhoZ) * length(vec_Vareps)
  mat_avgSim_eb_Normal <- matrix(ncol = 7,
                                 nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(mat_avgSim_eb_Normal) <- c("rhoZ", "Var_eps",
                                      paste0("theta", 1:p_Z),
                                      "ATT")
  mat_avgSim_ebme_Normal <- mat_avgSim_eb_Normal
  mat_avgSim_bc_Normal <- mat_avgSim_eb_Normal
  mat_avgSim_hcc_Normal <- mat_avgSim_eb_Normal
  mat_avgSim_hyj_Normal <- mat_avgSim_eb_Normal

  ### mat_BiasSdMSE:
  #### 17 cols: rhoZ, Var(eps), bias(theta1, theta2, theta3, theta4, ATT), sd(theta1, theta2, theta3, theta4, ATT), MSE(theta1, theta2, theta3, theta4, ATT)
  #### nrow = length(vec_rhoZ) * length(vec_Vareps)
  mat_BiasSdMSE_eb_Normal <- matrix(ncol = 17,
                                    nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(mat_BiasSdMSE_eb_Normal) <- c("rhoZ", "Var_eps",
                                         paste0("BIAS_", c(paste0("theta", 1:4), "ATT")),
                                         paste0("SD_", c(paste0("theta", 1:4), "ATT")),
                                         paste0("MSE_", c(paste0("theta", 1:4), "ATT")))
  mat_BiasSdMSE_ebme_Normal <- mat_BiasSdMSE_eb_Normal
  mat_BiasSdMSE_bc_Normal <- mat_BiasSdMSE_eb_Normal
  mat_BiasSdMSE_hcc_Normal <- mat_BiasSdMSE_eb_Normal
  mat_BiasSdMSE_hyj_Normal <- mat_BiasSdMSE_eb_Normal
  # mat_BiasSdMSE_oracle_Normal <- mat_BiasSdMSE_eb_Normal

  ## Iteration numbers and indicator pool
  ### EBME-BFGS
  iterbfgsind.r_Normal <- rep(TRUE, N_Sim)
  ### BC - approximate bias corrector
  iterbcind.r_Normal <- iterbfgsind.r_Normal
  ### HCC
  iterhccind.r_Normal <- iterbfgsind.r_Normal
  ### HYJ
  iterhyjind.r_Normal <- iterbfgsind.r_Normal

  ## Rate of model convergence: Calculate the rate of convergence in each combination of rho_Z and var_eps
  ### 3 cols: rhoZ, Var(eps), average value of the rate of convergence
  ### nrow <- length(vec_rhoZ) * length(vec_Vareps)
  Rc_ebme_Normal <- matrix(data = 0,
                           ncol = 3,
                           nrow = length(vec_rhoZ) * length(vec_Vareps))
  colnames(Rc_ebme_Normal) <- c("rhoZ", "Var_eps", "Rc")
  Rc_bc_Normal <- Rc_ebme_Normal
  Rc_hcc_Normal <- Rc_ebme_Normal
  Rc_hyj_Normal <- Rc_ebme_Normal

  # matrix markdown the ASMD
  mat_ASMD_CEB_Normal <- NULL
  mat_ASMD_naiveEB_Normal <- NULL
  mat_ASMD_BCEB_Normal <- NULL
  mat_ASMD_HL_Normal <- NULL
  mat_ASMD_HW_Normal <- NULL
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
      ## Normal distributed measurement errors
      mat_SingleSim_eb_Normal[i, "rhoZ"] <- rhoZ
      mat_SingleSim_ebme_Normal[i, "rhoZ"] <- rhoZ
      mat_SingleSim_bc_Normal[i, "rhoZ"] <- rhoZ
      mat_SingleSim_hcc_Normal[i, "rhoZ"] <- rhoZ
      mat_SingleSim_hyj_Normal[i, "rhoZ"] <- rhoZ
    }

      # Mark the Var_eps
    {
      ## Normal distributed measurement errors
      mat_SingleSim_eb_Normal[i, "Var_eps"] <- var_eps
      mat_SingleSim_ebme_Normal[i, "Var_eps"] <- var_eps
      mat_SingleSim_bc_Normal[i, "Var_eps"] <- var_eps
      mat_SingleSim_hcc_Normal[i, "Var_eps"] <- var_eps
      mat_SingleSim_hyj_Normal[i, "Var_eps"] <- var_eps
    }

      # True covariates
    {
      Z <- mvrnorm(n = N, mu = c(Mean_X, Mean_U), Sigma = Sigma_Z)
      trueX <- Z[, 1:p_X]
      U <- Z[, (p_X + 1):p_Z]
    }

      # Construct replicated covariate dataset
    {
      # Normal distributed measurement errors for replicated dataset
      for (k in 1:r) {
        measure.error_Normal[[k]] <- mvrnorm(n = N, mu = rep(0, p_X), Sigma = Sigma_e)

        X_Normal[[k]] <- trueX + measure.error_Normal[[k]]
        Ze_Normal[[k]] <- data.frame(cbind(X_Normal[[k]], U))
        colnames(Ze_Normal[[k]]) <- c("X1", "X2", "U3", "U4")
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
      ## Normal distributed measurement errors
      Ze1_bar_al_Normal <- rep(0, p_Z)
      for (k in 1:r) {
        Ze0_Normal[[k]] <- as.matrix(Ze_Normal[[k]][Index_col,])
        Ze1_Normal[[k]] <- as.matrix(Ze_Normal[[k]][Index_trt,])
        Ze1_bar_Normal[[k]] <- colMeans(Ze1_Normal[[k]])
        Ze1_bar_al_Normal <- Ze1_bar_al_Normal + Ze1_bar_Normal[[k]]
        Ze_Normal[[k]] <- as.matrix(Ze_Normal[[k]])
      }
      Ze1_bar_al_Normal <- Ze1_bar_al_Normal / r
    }

      # Estimate the variance-covariance matrix $\Sigma$ of measurement error
    {
      meanZe_al_Normal <- Reduce("+", Ze_Normal) / r
      erhat_Normal <- matrix(rep(0, ncol(Z)^2), ncol = ncol(Z), nrow = ncol(Z))
      for (kk in 1:r) {
        erhat_Normal <- erhat_Normal +
          t(Ze_Normal[[kk]] - meanZe_al_Normal) %*% (Ze_Normal[[kk]] - meanZe_al_Normal)
      }
      Sigmahat_al_Normal <- erhat_Normal / ((r - 1) * N)
      # round(Sigmahat_al_Normal, 3)
    }

      # ZQY-EB: native Entropy Balancing
    {
      opt.out.zqy.r_Normal <- optim(par = rep(0, (ncol(Z))),
                                    fn = objective.zqy.r_Normal,
                                    gr = gradient.zqy.r_Normal,
                                    method = "BFGS",
                                    control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200))
      theta0_Normal <- opt.out.zqy.r_Normal$par
      # round(theta0_Normal, 3)
    }

      # EBME - BFGS
    {
      output.CEB.solveObj_Normal <- tryCatch(optim(par = rep(0, p_Z),
                                                   fn = objective.ebme.r_Normal,
                                                   gr = gradient.ebme.r_Normal,
                                                   method = "BFGS",
                                                   control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200)),
                                             error = function(e) {
                                               print("Fail to solve objective of CEB, We directly use naive EB's results for subsequent comparison.")
                                               return(FALSE)
                                             })
      if (class(output.CEB.solveObj_Normal) != "list" |
        is.na(output.CEB.solveObj_Normal$value) |
        is.nan(output.CEB.solveObj_Normal$value)) {
        output.CEB.solveObj_Normal <- opt.out.zqy.r_Normal
      }
      theta_CEB_solveObj <- output.CEB.solveObj_Normal$par
      ObjValue_CEB_solveObj <- output.CEB.solveObj_Normal$value

      # theta_CEB - solve the gradient funtion
      output.CEB.solveGradient_Normal <- tryCatch(multiroot(f = gradient.ebme.r_Normal,
                                                            start = rep(0, p_Z),
                                                            rtol = .Machine$double.eps, maxiter = 200),
                                                  error = function(e) {
                                                    print("Fail to solve gradient of CEB, We directly use naive EB's results for subsequent comparison.")
                                                    return(FALSE)
                                                  })
      if (class(output.CEB.solveGradient_Normal) != "list") {
        ObjValue_CEB_solveGradient <- 9999
      }else {
        theta_CEB_solveGradient <- output.CEB.solveGradient_Normal$root
        ObjValue_CEB_solveGradient <- objective.ebme.r_Normal(theta = theta_CEB_solveGradient)
        if (is.nan(ObjValue_CEB_solveGradient) | is.na(ObjValue_CEB_solveGradient)) { ObjValue_CEB_solveGradient <- 9999 }
      }

      if (ObjValue_CEB_solveObj <= ObjValue_CEB_solveGradient) {
        theta_CEB_Normal <- theta_CEB_solveObj
        iterbfgsind.r_Normal[i] <- (max(abs(gradient.ebme.r_Normal(theta_CEB_solveObj))) < 1e-04)
        ## double check
        if (is.na(iterbfgsind.r_Normal[i]) | is.nan(iterbfgsind.r_Normal[i])) {
          iterbfgsind.r_Normal[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }else {
        theta_CEB_Normal <- theta_CEB_solveGradient
        iterbfgsind.r_Normal[i] <- (max(abs(output.CEB.solveGradient_Normal$f.root)) < sqrt(.Machine$double.eps)) & (output.CEB.solveGradient_Normal$estim.precis != "NaN")
        ## double check
        if (is.na(iterbfgsind.r_Normal[i]) | is.nan(iterbfgsind.r_Normal[i])) {
          iterbfgsind.r_Normal[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }
    }

      # BCEB
    {
      ## Normal distributed measurement errors
      theta_bc_Normal <- drop(ginv(hessian.zqy.r_Normal(theta0_Normal) - Sigma) %*%
                                hessian.zqy.r_Normal(theta0_Normal) %*%
                                theta0_Normal)
      iterbcind.r_Normal[i] <- !(NaN %in% theta_bc_Normal | (sum(is.infinite(theta_bc_Normal)) > 0))
    }

      # HCC(2012)
    {
      ## Normal distributed measurement errors
      opt.out.hcc.r_Normal <- multiroot(f = gradient.hcc.r_Normal, start = rep(0, p_Z),
                                        rtol = .Machine$double.eps, maxiter = 200)
      iterhccind.r_Normal[i] <- (max(gradient.hcc.r_Normal(opt.out.hcc.r_Normal$root)) < sqrt(.Machine$double.eps)) & (opt.out.hcc.r_Normal$estim.precis != "NaN")
    }

      # HYJ(2012)
    {
      ## Normal distributed measurement errors
      opt.out.hyj.r_Normal <- multiroot(f = gradient.hyj.r_Normal, start = rep(0, p_Z),
                                        rtol = .Machine$double.eps, maxiter = 200)
      iterhyjind.r_Normal[i] <- (max(gradient.hyj.r_Normal(opt.out.hyj.r_Normal$root)) < sqrt(.Machine$double.eps)) & (opt.out.hyj.r_Normal$estim.precis != "NaN")
      ### double check
      if (is.na(iterhyjind.r_Normal[i])) {
        iterhyjind.r_Normal[i] <- FALSE
        print("NA clear!")
      }
    }

      # ATT estimation with \hat{\theta}
    {
      # ZQY-EB: naive Entropy Balancing
    {
      weight_zqy_Normal <- exp(Ze0_Normal[[1]] %*% theta0_Normal)
      weight_zqy_Normal <- weight_zqy_Normal / sum(weight_zqy_Normal)
      ATT_zqy_Normal <- mean(Y1) - t(weight_zqy_Normal) %*% Y0
      mat_SingleSim_eb_Normal[i, "ATT"] <- mat_SingleSim_eb_Normal[i, "ATT"] + ATT_zqy_Normal
    }

      # CEB / EBME - bfgs C2: the proposed method
    {
      weight_bfgs_Normal <- exp(Ze0_Normal[[1]] %*% theta_CEB_Normal)
      weight_bfgs_Normal <- weight_bfgs_Normal / sum(weight_bfgs_Normal)
      ATT_bfgs_Normal <- mean(Y1) - t(weight_bfgs_Normal) %*% Y0
      mat_SingleSim_ebme_Normal[i, "ATT"] <- mat_SingleSim_ebme_Normal[i, "ATT"] + ATT_bfgs_Normal
    }

      # BCEB
    {
      weight_bc_Normal <- exp(Ze0_Normal[[1]] %*% theta_bc_Normal)
      weight_bc_Normal <- weight_bc_Normal / sum(weight_bc_Normal)
      ATT_bc_Normal <- mean(Y1) - t(weight_bc_Normal) %*% Y0
      mat_SingleSim_bc_Normal[i, "ATT"] <- mat_SingleSim_bc_Normal[i, "ATT"] + ATT_bc_Normal
    }

      # HCC
    {
      par_hcc_Normal <- as.matrix(opt.out.hcc.r_Normal$root, ncol = 1)
      weight_hcc_Normal <- matrix(rep(0, nrow(weight_bc_Normal)))
      for (k in 1:r) {
        weight_hcc_Normal <- weight_hcc_Normal + exp(Ze0_Normal[[k]] %*% par_hcc_Normal)
      }
      weight_hcc_Normal <- weight_hcc_Normal / sum(weight_hcc_Normal)
      ATT_hcc_Normal <- mean(Y1) - t(weight_hcc_Normal) %*% Y0
      mat_SingleSim_hcc_Normal[i, "ATT"] <- mat_SingleSim_hcc_Normal[i, "ATT"] + ATT_hcc_Normal
    }

      # HYJ
    {
      par_hyj_Normal <- as.matrix(opt.out.hyj.r_Normal$root, ncol = 1)
      weight_hyj_Normal <- matrix(rep(0, nrow(weight_bc_Normal)))
      for (k in 1:r) {
        weight_hyj_Normal <- weight_hyj_Normal + exp(Ze0_Normal[[k]] %*% par_hyj_Normal)
      }
      weight_hyj_Normal <- weight_hyj_Normal / sum(weight_hyj_Normal)
      ATT_hyj_Normal <- mean(Y1) - t(weight_hyj_Normal) %*% Y0
      ### double check
      if (is.na(ATT_hyj_Normal)) {
        iterhyjind.r_Normal[i] <- FALSE
        print("ATT of HW - NA clear!")
      }
      mat_SingleSim_hyj_Normal[i, "ATT"] <- mat_SingleSim_hyj_Normal[i, "ATT"] + ATT_hyj_Normal
    }

    }

      # theta coefficient
    {
      mat_SingleSim_eb_Normal[i, paste0("theta", 1:p_Z)] <- theta0_Normal
      mat_SingleSim_ebme_Normal[i, paste0("theta", 1:p_Z)] <- theta_CEB_Normal
      mat_SingleSim_bc_Normal[i, paste0("theta", 1:p_Z)] <- theta_bc_Normal
      mat_SingleSim_hcc_Normal[i, paste0("theta", 1:p_Z)] <- opt.out.hcc.r_Normal$root
      mat_SingleSim_hyj_Normal[i, paste0("theta", 1:p_Z)] <- opt.out.hyj.r_Normal$root
    }

      # ASMD calculation ----
    {
      # mat_SingleSim_ASMD_CEB
      mat_SingleSim_ASMD_naiveEB[i,] <- unlist(ASMDCal(weight_zqy_Normal))
      mat_SingleSim_ASMD_CEB[i,] <- unlist(ASMDCal(weight_bfgs_Normal))
      mat_SingleSim_ASMD_BCEB[i,] <- unlist(ASMDCal(weight_bc_Normal))
      mat_SingleSim_ASMD_HL[i,] <- unlist(ASMDCal(weight_hcc_Normal))
      mat_SingleSim_ASMD_HW[i,] <- unlist(ASMDCal(weight_hyj_Normal))
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
    Rc_ebme_Normal[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterbfgsind.r_Normal) / N_Sim)
    ## BC
    Rc_bc_Normal[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterbcind.r_Normal) / N_Sim)
    ## HCC
    Rc_hcc_Normal[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterhccind.r_Normal) / N_Sim)
    ## HYJ
    Rc_hyj_Normal[ith_vareps_rhoZ,] <- c(rhoZ, var_eps, sum(iterhyjind.r_Normal) / N_Sim)
  }

    # mat_avgSim赋值 - Convergent results
  {
    mat_avgSim_eb_Normal[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_eb_Normal)
    mat_avgSim_ebme_Normal[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_ebme_Normal[iterbfgsind.r_Normal,])
    mat_avgSim_bc_Normal[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_bc_Normal[iterbcind.r_Normal,])
    mat_avgSim_hcc_Normal[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_hcc_Normal[iterhccind.r_Normal,])
    mat_avgSim_hyj_Normal[ith_vareps_rhoZ,] <- colMeans(mat_SingleSim_hyj_Normal[iterhyjind.r_Normal,])
  }

    # mat_BiasSDMSE - Convergent results
    trueVal <- c(True.theta[-1], ATE)
  {
    ## rhoZ
    mat_BiasSdMSE_eb_Normal[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_ebme_Normal[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_bc_Normal[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_hcc_Normal[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    mat_BiasSdMSE_hyj_Normal[ith_vareps_rhoZ, "rhoZ"] <- rhoZ
    ## Var_eps
    mat_BiasSdMSE_eb_Normal[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_ebme_Normal[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_bc_Normal[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_hcc_Normal[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    mat_BiasSdMSE_hyj_Normal[ith_vareps_rhoZ, "Var_eps"] <- var_eps
    ## bias
    ColnamesinNeed_avgmat <- c(paste0("theta", 1:p_Z), "ATT")
    ColnamesinNeed_Bias <- paste0("BIAS_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Normal[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_eb_Normal[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_ebme_Normal[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_ebme_Normal[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_bc_Normal[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_bc_Normal[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_hcc_Normal[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_hcc_Normal[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    mat_BiasSdMSE_hyj_Normal[ith_vareps_rhoZ, ColnamesinNeed_Bias] <- mat_avgSim_hyj_Normal[ith_vareps_rhoZ, ColnamesinNeed_avgmat] - trueVal
    ## sd
    ColnamesinNeed_SD <- paste0("SD_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Normal[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_eb_Normal[, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_ebme_Normal[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_ebme_Normal[iterbfgsind.r_Normal, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_bc_Normal[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_bc_Normal[iterbcind.r_Normal, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_hcc_Normal[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_hcc_Normal[iterhccind.r_Normal, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    mat_BiasSdMSE_hyj_Normal[ith_vareps_rhoZ, ColnamesinNeed_SD] <- apply(mat_SingleSim_hyj_Normal[iterhyjind.r_Normal, ColnamesinNeed_avgmat], MARGIN = 2, FUN = sd)
    ## MSE
    ColnamesinNeed_MSE <- paste0("MSE_", c(paste0("theta", 1:p_Z), "ATT"))
    mat_BiasSdMSE_eb_Normal[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_eb_Normal[, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_ebme_Normal[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_ebme_Normal[iterbfgsind.r_Normal, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_bc_Normal[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_bc_Normal[iterbcind.r_Normal, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_hcc_Normal[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_hcc_Normal[iterhccind.r_Normal, ColnamesinNeed_avgmat])
    mat_BiasSdMSE_hyj_Normal[ith_vareps_rhoZ, ColnamesinNeed_MSE] <- MSE(mat_SingleSim_hyj_Normal[iterhyjind.r_Normal, ColnamesinNeed_avgmat])
  }

    # ASMD
  {
    # mat_ASMD_CEB_Normal
    mat_ASMD_naiveEB_Normal <- rbind(mat_ASMD_naiveEB_Normal,
                                     colMeans(mat_SingleSim_ASMD_naiveEB, na.rm = TRUE))
    mat_ASMD_CEB_Normal <- rbind(mat_ASMD_CEB_Normal,
                                 colMeans(mat_SingleSim_ASMD_CEB, na.rm = TRUE))
    mat_ASMD_BCEB_Normal <- rbind(mat_ASMD_BCEB_Normal,
                                  colMeans(mat_SingleSim_ASMD_BCEB, na.rm = TRUE))
    mat_ASMD_HL_Normal <- rbind(mat_ASMD_HL_Normal,
                                colMeans(mat_SingleSim_ASMD_HL, na.rm = TRUE))
    mat_ASMD_HW_Normal <- rbind(mat_ASMD_HW_Normal,
                                colMeans(mat_SingleSim_ASMD_HW, na.rm = TRUE))
    # round(mat_ASMD_CEB_Normal, 3)
  }

    Count_Row <- 0
  {
  {
    # mat_SingleSim:
    ## 7 cols: rhoZ, Var(eps), theta1, theta2, theta3, theta4, ATT
    ## nrow = N_Sim
    mat_SingleSim_eb_Normal <- matrix(data = 0,
                                      ncol = 7, nrow = N_Sim)
    colnames(mat_SingleSim_eb_Normal) <- c("rhoZ", "Var_eps",
                                           paste0("theta", 1:p_Z), "ATT")
    mat_SingleSim_ebme_Normal <- mat_SingleSim_eb_Normal
    mat_SingleSim_bc_Normal <- mat_SingleSim_eb_Normal
    mat_SingleSim_hcc_Normal <- mat_SingleSim_eb_Normal
    mat_SingleSim_hyj_Normal <- mat_SingleSim_eb_Normal

    ## EBME-BFGS
    iterbfgsind.r_Normal <- rep(TRUE, N_Sim)
    ## BC - approximate bias corrector
    iterbcind.r_Normal <- iterbfgsind.r_Normal
    ## HCC
    iterhccind.r_Normal <- iterbfgsind.r_Normal
    ## HYJ
    iterhyjind.r_Normal <- iterbfgsind.r_Normal

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
save.image(file = "CompareAllMethods_NormalEps_Original.RData")

{
  # Table output----
{
  Rc_ebme_Normal <- cbind(Rc_ebme_Normal, Method = "CEB")
  Rc_bc_Normal <- cbind(Rc_bc_Normal, Method = "BCEB")
  Rc_hcc_Normal <- cbind(Rc_hcc_Normal, Method = "HL")
  Rc_hyj_Normal <- cbind(Rc_hyj_Normal, Method = "HW")
  mat_Rc_Normal <- rbind(Rc_ebme_Normal,
                         Rc_bc_Normal,
                         Rc_hcc_Normal,
                         Rc_hyj_Normal)
  csvName_mat_Rc_Normal <- paste0("mat_NormalEps_ConvergenceRate_N", N,
                                  "_Sim", N_Sim, ".csv")
  write.table(mat_Rc_Normal, file = csvName_mat_Rc_Normal,
              sep = ",", row.names = F)

  ## mat_avgSim
  mat_avgSim_eb_Normal <- cbind(mat_avgSim_eb_Normal, Method = "Naive")
  mat_avgSim_ebme_Normal <- cbind(mat_avgSim_ebme_Normal, Method = "CEB")
  mat_avgSim_bc_Normal <- cbind(mat_avgSim_bc_Normal, Method = "BCEB")
  mat_avgSim_hcc_Normal <- cbind(mat_avgSim_hcc_Normal, Method = "HL")
  mat_avgSim_hyj_Normal <- cbind(mat_avgSim_hyj_Normal, Method = "HW")
  mat_avgSim_Normal <- rbind(mat_avgSim_eb_Normal,
                             mat_avgSim_ebme_Normal,
                             mat_avgSim_bc_Normal,
                             mat_avgSim_hcc_Normal,
                             mat_avgSim_hyj_Normal)
  csvName_mat_avgSim <- paste0("mat_NormalEps_avgSim_theta_ATT_N", N, "_Sim", N_Sim,
                               ".csv")
  write.table(mat_avgSim_Normal, file = csvName_mat_avgSim,
              sep = ",", row.names = F)

  ## mat_BiasSdMSE
  mat_BiasSdMSE_eb_Normal <- cbind(mat_BiasSdMSE_eb_Normal, Method = "Naive")
  mat_BiasSdMSE_ebme_Normal <- cbind(mat_BiasSdMSE_ebme_Normal, Method = "CEB")
  mat_BiasSdMSE_bc_Normal <- cbind(mat_BiasSdMSE_bc_Normal, Method = "BCEB")
  mat_BiasSdMSE_hcc_Normal <- cbind(mat_BiasSdMSE_hcc_Normal, Method = "HL")
  mat_BiasSdMSE_hyj_Normal <- cbind(mat_BiasSdMSE_hyj_Normal, Method = "HW")
  mat_BiasSdMSE_Normal <- rbind(mat_BiasSdMSE_eb_Normal,
                                mat_BiasSdMSE_ebme_Normal,
                                mat_BiasSdMSE_bc_Normal,
                                mat_BiasSdMSE_hcc_Normal,
                                mat_BiasSdMSE_hyj_Normal)
  csvName_mat_BiasSdMSE <- paste0("mat_NormalEps_BiasSdMSE_theta_ATT_N", N, "_Sim", N_Sim,
                                  ".csv")
  write.table(mat_BiasSdMSE_Normal, file = csvName_mat_BiasSdMSE,
              sep = ",", row.names = F)
}

  # Plots----
{
  # Plot of convergenc rate of EBME, BC, HCC, HYJ
  Rcdf_Normal <- as.data.frame(mat_Rc_Normal)
  # Rcdf_Normal$Method <- factor(Rcdf_Normal$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  Rcdf_Normal$Method <- factor(Rcdf_Normal$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  p_Rc_NormalEps <- ggplot(Rcdf_Normal) +
    geom_line(aes(x = Var_eps, y = Rc, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = Rc, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Rate of Convergence", title = "Measurement Error with a Normal Distribution") +
    theme_bw()
  p_Rc_NormalEps
  ggsave(filename = "p_NormalEps_Rc.pdf", p_Rc_NormalEps,
         width = 14, height = 13,
         units = "cm")

  ## Seperate bias, sd, MSE dataframe
  mat_bias_Normal <- mat_BiasSdMSE_Normal[, c("rhoZ", "Var_eps", paste0("BIAS_", c(paste0("theta", 1:4), "ATT")), "Method")]
  mat_sd_Normal <- mat_BiasSdMSE_Normal[, c("rhoZ", "Var_eps", paste0("SD_", c(paste0("theta", 1:4), "ATT")), "Method")]
  mat_MSE_Normal <- mat_BiasSdMSE_Normal[, c("rhoZ", "Var_eps", paste0("MSE_", c(paste0("theta", 1:4), "ATT")), "Method")]
  ## bias
  biasdf_Normal <- melt(as.data.frame(mat_bias_Normal),
                        id.vars = c("rhoZ", "Var_eps", "Method"))
  biasdf_Normal$value <- as.numeric(biasdf_Normal$value)
  which(is.na(biasdf_Normal$value))
  # biasdf_Normal$Method <- factor(biasdf_Normal$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  biasdf_Normal$Method <- factor(biasdf_Normal$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### bias_theta
  biasdf_theta_Normal <- biasdf_Normal[biasdf_Normal$variable != "BIAS_ATT",]
  p_AllMethods_bias_theta_NormalEps <- ggplot(biasdf_theta_Normal) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Bias", title = "True Value: theta1 = -1, theta2 = 0.5, theta2 = -0.25, theta4 = -0.1") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_bias_theta_NormalEps
  ggsave(filename = "p_NormalEps_AllMethods_bias_theta.pdf", p_AllMethods_bias_theta_NormalEps,
         width = 26, height = 26,
         units = "cm")
  ### bias_ATT
  biasdf_ATT_Normal <- biasdf_Normal[biasdf_Normal$variable == "BIAS_ATT",]
  p_AllMethods_bias_ATT_NormalEps <- ggplot(biasdf_ATT_Normal) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Bias") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_bias_ATT_NormalEps
  ggsave(filename = "p_NormalEps_AllMethods_bias_ATT.pdf", p_AllMethods_bias_ATT_NormalEps,
         width = 26, height = 10,
         units = "cm")
  ## sd
  sddf_Normal <- melt(as.data.frame(mat_sd_Normal),
                      id.vars = c("rhoZ", "Var_eps", "Method"))
  sddf_Normal$value <- as.numeric(sddf_Normal$value)
  which(is.na(sddf_Normal$value))
  # sddf_Normal$Method <- factor(sddf_Normal$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  sddf_Normal$Method <- factor(sddf_Normal$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### sd_theta
  sddf_theta_Normal <- sddf_Normal[sddf_Normal$variable != "SD_ATT",]
  p_AllMethods_sd_theta_NormalEps <- ggplot(sddf_theta_Normal) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Standard Deviation", title = "SD of Coefficients with Normal Distributed Measurement Error") +
    theme_bw()
  p_AllMethods_sd_theta_NormalEps
  ggsave(filename = "p_NormalEps_AllMethods_sd_theta.pdf", p_AllMethods_sd_theta_NormalEps,
         width = 26, height = 26,
         units = "cm")
  ### sd_ATT
  sddf_ATT_Normal <- sddf_Normal[sddf_Normal$variable == "SD_ATT",]
  p_AllMethods_sd_ATT_NormalEps <- ggplot(sddf_ATT_Normal) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    scale_y_continuous(limits = c(0, max(sddf_ATT_Normal$value)), breaks = seq(from = 0, to = ceiling(max(sddf_ATT_Normal$value)), by = 0.5)) +
    labs(x = "Variance of Measurement Error", y = "Standard Deviation") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_sd_ATT_NormalEps
  ggsave(filename = "p_NormalEps_AllMethods_sd_ATT.pdf", p_AllMethods_sd_ATT_NormalEps,
         width = 26, height = 10,
         units = "cm")
  ## MSE
  MSEdf_Normal <- melt(as.data.frame(mat_MSE_Normal),
                       id.vars = c("rhoZ", "Var_eps", "Method"))
  MSEdf_Normal$value <- as.numeric(MSEdf_Normal$value)
  which(is.na(MSEdf_Normal$value))
  # MSEdf_Normal$Method <- factor(MSEdf_Normal$Method, levels = c("EB", "EBME", "BC", "HCC", "HYJ"))
  MSEdf_Normal$Method <- factor(MSEdf_Normal$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ### MSE_theta
  MSEdf_theta_Normal <- MSEdf_Normal[MSEdf_Normal$variable != "MSE_ATT",]
  p_AllMethods_MSE_theta_NormalEps <- ggplot(MSEdf_theta_Normal) +
    facet_grid(rows = vars(variable), cols = vars(rhoZ)) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Mean Squared Error", title = "MSE of Coefficients with Normal Distributed Measurement Error") +
    theme_bw()
  p_AllMethods_MSE_theta_NormalEps
  ggsave(filename = "p_NormalEps_AllMethods_MSE_theta.pdf", p_AllMethods_MSE_theta_NormalEps,
         width = 26, height = 26,
         units = "cm")
  ### MSE_ATT
  MSEdf_ATT_Normal <- MSEdf_Normal[MSEdf_Normal$variable == "MSE_ATT",]
  p_AllMethods_MSE_ATT_NormalEps <- ggplot(MSEdf_ATT_Normal) +
    geom_line(aes(x = Var_eps, y = value, color = Method, group = Method)) +
    geom_point(aes(x = Var_eps, y = value, color = Method, shape = Method)) +
    labs(x = "Variance of Measurement Error", y = "Mean Squared Error") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    theme_bw()
  p_AllMethods_MSE_ATT_NormalEps
  ggsave(filename = "p_NormalEps_AllMethods_MSE_ATT.pdf", p_AllMethods_MSE_ATT_NormalEps,
         width = 26, height = 10,
         units = "cm")

  Figure_NormalEps_BiasSDMSE_ATT <- ggarrange(p_AllMethods_bias_ATT_NormalEps,
                                              p_AllMethods_sd_ATT_NormalEps,
                                              p_AllMethods_MSE_ATT_NormalEps,
                                              nrow = 1, ncol = 3,
                                              widths = c(8, 8, 8, 8), labels = c('(a)', '(b)', '(c)'),
                                              vjust = 1.1, hjust = 0,
                                              common.legend = T, legend = 'right',
                                              font.label = list(size = 23, face = 'plain'))
  Figure_NormalEps_BiasSDMSE_ATT
}
}

# Save .Rdata
save.image(file = "CompareAllMethods_NormalEps_avgEdited.RData")

# Figure Plot: Combine the plot of normal distributed and Normal distribution----
{
  # Plot of bias+SD+MSE of ATT - Normal
  Figure_NormalEps_BiasSDMSE_ATT <- ggarrange(p_AllMethods_bias_ATT_NormalEps,
                                              p_AllMethods_sd_ATT_NormalEps,
                                              p_AllMethods_MSE_ATT_NormalEps,
                                              # nrow = 3, ncol = 3,
                                              nrow = 1, ncol = 3,
                                              widths = c(8, 8, 8, 8, 8),
                                              # heights = c(8, 8, 8),
                                              heights = 8,
                                              common.legend = T, legend = 'top',
                                              labels = c('Case I', '', ''),
                                              vjust = 0, hjust = -0.6,
                                              font.label = list(size = 18, face = 'bold'))
  Figure_NormalEps_BiasSDMSE_ATT

  ggsave(filename = "p_NormalEps_BiasSDMSE_ATT.pdf",
         Figure_NormalEps_BiasSDMSE_ATT,
         width = 33, height = 12,
         units = "cm")
}

# -------------------------------------- #
mat_ASMD_CEB_Normal <- cbind(var_eps = vec_Vareps, mat_ASMD_CEB_Normal)
colnames(mat_ASMD_CEB_Normal)[-1] <- paste0(colnames(Ze_Normal[[1]]), "ASMD")

mat_ASMD_naiveEB_Normal <- cbind(var_eps = vec_Vareps, mat_ASMD_naiveEB_Normal)
colnames(mat_ASMD_naiveEB_Normal)[-1] <- paste0(colnames(Ze_Normal[[1]]), "ASMD")

mat_ASMD_BCEB_Normal <- cbind(var_eps = vec_Vareps, mat_ASMD_BCEB_Normal)
colnames(mat_ASMD_BCEB_Normal)[-1] <- paste0(colnames(Ze_Normal[[1]]), "ASMD")

mat_ASMD_HL_Normal <- cbind(var_eps = vec_Vareps, mat_ASMD_HL_Normal)
colnames(mat_ASMD_HL_Normal)[-1] <- paste0(colnames(Ze_Normal[[1]]), "ASMD")

mat_ASMD_HW_Normal <- cbind(var_eps = vec_Vareps, mat_ASMD_HW_Normal)
colnames(mat_ASMD_HW_Normal)[-1] <- paste0(colnames(Ze_Normal[[1]]), "ASMD")


round(mat_ASMD_naiveEB_Normal, 3)
round(mat_ASMD_CEB_Normal, 3)
round(mat_ASMD_BCEB_Normal, 3)
round(mat_ASMD_HL_Normal, 3)
round(mat_ASMD_HW_Normal, 3)

# Export the data of ASMD into the same Excel with different sheets
library(openxlsx)

list_of_ASMD_NormalME <- list("naiveEB" = round(mat_ASMD_naiveEB_Normal, 3),
                              "CEB" = round(mat_ASMD_CEB_Normal, 3),
                              "BCEB" = round(mat_ASMD_BCEB_Normal, 3),
                              "HL" = round(mat_ASMD_HL_Normal, 3),
                              "HW" = round(mat_ASMD_HW_Normal, 3))

openxlsx::write.xlsx(x = list_of_ASMD_NormalME,
                     file = "list_of_ASMD_NormalME.xlsx",
                     row.names = FALSE)

mat_ASMD_naiveEB_Normal_1 <- cbind(Method = "naiveEB", t(round(mat_ASMD_naiveEB_Normal, 3)))
mat_ASMD_CEB_Normal_1 <- cbind(Method = "CEB", t(round(mat_ASMD_CEB_Normal, 3)))
mat_ASMD_BCEB_Normal_1 <- cbind(Method = "BCEB", t(round(mat_ASMD_BCEB_Normal, 3)))
mat_ASMD_HL_Normal_1 <- cbind(Method = "HL", t(round(mat_ASMD_HL_Normal, 3)))
mat_ASMD_HW_Normal_1 <- cbind(Method = "HW", t(round(mat_ASMD_HW_Normal, 3)))

# 一定要转成dataframe才能成功输出
mat_ASMD_AllMethods <- as.data.frame(rbind(mat_ASMD_naiveEB_Normal_1,
                                           mat_ASMD_CEB_Normal_1,
                                           mat_ASMD_BCEB_Normal_1,
                                           mat_ASMD_HL_Normal_1,
                                           mat_ASMD_HW_Normal_1))

openxlsx::write.xlsx(x = mat_ASMD_AllMethods,
                     file = "mat_ASMD_AllMethods_NormalME.xlsx",
                     row.names = TRUE)