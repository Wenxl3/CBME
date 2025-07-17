# -*- coding: utf-8 -*-
# ------------------------------------------
# Project     : EBME_R_202304242131
# Title       : ASMD Analysis of CEB
# Objective   : Analyze the performance of CEB with normal measurement error.
# Created on  : 2023/05/09 10:52
# ------------------------------------------

{
  rm(list = ls())
  getwd()
  # To suppress warnings in the global settings
  options(warn = -1)
}

# Importing Self-Defined Function----
## Objective and gradient function
{

  # ======================================== #
  #           Method of naive EB             #
  # ======================================== #

  # EB: Entropy Balancing (Jens Hainmueller)----
  objective.zqy.r_Normal <- function(theta) {
    f <- 0
    f <- f +
      log(sum(drop((1 - W) * exp(Ze_Normal[[1]] %*% theta)))) +
      sum(-Ze1_bar_Normal[[1]] * theta)
    return(f)
  }

  gradient.zqy.r_Normal <- function(theta) {
    f <- 0
    w <- drop((1 - W) * exp(Ze_Normal[[1]] %*% theta))
    f <- f + drop(w %*% (Ze_Normal[[1]]) / sum(w) - Ze1_bar_Normal[[1]])
    return(f)
  }

  hessian.zqy.r_Normal <- function(theta) {
    f <- 0
    w <- drop((1 - W) * exp(Ze_Normal[[1]] %*% theta))
    a <- matrix(rep(0, p_Z^2), nrow = p_Z, ncol = p_Z)
    for (j in 1:N) {
      a <- a + (Ze_Normal[[1]][j,]) %*% t(Ze_Normal[[1]][j,]) * w[j]
    }
    w2 <- a * sum(w) - t(w %*% Ze_Normal[[1]]) %*% (w %*% Ze_Normal[[1]])
    f <- f + drop(w2 / (sum(w))^2)
    return(f)
  }

  # ======================================== #
  #              Method of CEB               #
  # ======================================== #

  # CEB / EBME: our proposed method
  objective.ebme.r_Normal <- function(theta) {
    f <- drop(log((1 - W) %*% exp(Ze_Normal[[1]] %*% theta)) -
                Ze1_bar_Normal[[1]] %*% theta -
                t(theta) %*% Sigma %*% theta / 2)
    return(f)
  }

  gradient.ebme.r_Normal <- function(theta) {
    w <- drop((1 - W) * exp(Ze_Normal[[1]] %*% theta))
    f <- drop(w %*% Ze_Normal[[1]] / sum(w) - Ze1_bar_Normal[[1]]) -
      drop(Sigma %*% theta)
    return(f)
  }

  # CEB / EBME: our proposed method
  objective.ebme.r2_Normal <- function(theta) {
    f <- drop(log((1 - W) %*% exp(Ze_Normal[[2]] %*% theta)) -
                Ze1_bar_Normal[[2]] %*% theta -
                t(theta) %*% Sigma %*% theta / 2)
    return(f)
  }

  gradient.ebme.r2_Normal <- function(theta) {
    w <- drop((1 - W) * exp(Ze_Normal[[2]] %*% theta))
    f <- drop(w %*% Ze_Normal[[2]] / sum(w) - Ze1_bar_Normal[[2]]) -
      drop(Sigma %*% theta)
    return(f)
  }

  # ======================================== #
  #              Method of HL                #
  # ======================================== #

  # Chengcheng Hu & D. Y Lin (2004): replicate method 1----
  gradient.hcc.r_Normal <- function(theta) {
    eta0 <- 0
    eta1 <- 0
    for (u in 2:r) {
      for (s in 1:(u - 1)) {
        eta0 <- eta0 + sum(exp(drop((theta %*% t(Ze_Normal[[u]] - Ze_Normal[[s]])))) +
                             exp(drop((theta %*% t(Ze_Normal[[s]] - Ze_Normal[[u]])))))
        eta1 <- eta1 +
          (drop(t(Ze_Normal[[u]] - Ze_Normal[[s]]) %*% (exp(drop(theta %*% t(Ze_Normal[[u]] - Ze_Normal[[s]])))))) +
          (drop(t(Ze_Normal[[s]] - Ze_Normal[[u]]) %*% (exp(drop(theta %*% t(Ze_Normal[[s]] - Ze_Normal[[u]]))))))
      }
    }
    eta0 <- sqrt(eta0 / (r * (r - 1) * N))
    eta1 <- eta1 / (2 * N * r * (r - 1) * eta0)
    f0 <- 0
    f1 <- 0
    for (k in 1:r) {
      f0 <- f0 + sum(drop((1 - W) * exp(Ze_Normal[[k]] %*% theta)))
      f1 <- f1 + drop((t(Ze_Normal[[k]]) - eta1 / eta0) %*% drop((1 - W) * exp(Ze_Normal[[k]] %*% theta)))
    }
    f <- drop(f1 / f0 - Ze1_bar_al_Normal)
    # f <- sum(f^2)
    return(f)
  }

  # ======================================== #
  #              Method of HW                #
  # ======================================== #

  # Yijian Huang & C. Y. Wang (2000): replicate method 2----
  gradient.hyj.r_Normal <- function(theta) {
    index_Permutation <- permutations(n = r, r = 2, v = 1:r)
    N_Permutation <- nrow(index_Permutation)
    # RepDataset_Z <- list(repset1 = Ze_Normal[[1]], repset2 = Ze_Normal[[2]], repset3 = Ze_Normal[[3]])
    RepDataset_Z <- list(repset1 = Ze_Normal[[1]], repset2 = Ze_Normal[[2]])

    for (i in 1:N_Permutation) {
      # i <- 1

      # Permutaition Set Construction
      ## Set_Permutation_i: the ith permutaiton of replicated dataset
      ## Set_Permutation_i[[j]]: the jth covariate matrix in the ith permutaiton of replicated dataset(i=1,…,N_Permutation, N_j = 1,2)
      Set_Permutation_i <- RepDataset_Z[index_Permutation[i,]]

      # 将所有permutaion组合下 exp(theta * Z1)*Z2 的运算结果求和
      if (i == 1) {
        # mat_Permutation: 记录所有permutaion组合下 exp(theta * Z1)*Z2 的运算结果
        mat_Permutation <- drop(exp(Set_Permutation_i[[1]] %*% theta)) * Set_Permutation_i[[2]]
      }else {
        mat_Permutation <- mat_Permutation + drop(exp(Set_Permutation_i[[1]] %*% theta)) * Set_Permutation_i[[2]]
      }
    }
    avgmat_Permutation <- mat_Permutation / r / (r - 1)
    Numerator <- colSums((1 - W) * avgmat_Permutation)
    # 注意：simulation比较特殊的设定是每n_i都相同，n_i = r
    for (k in 1:r) {
      if (k == 1) {
        Denominator <- (1 - W) %*% exp(Ze_Normal[[k]] %*% theta)
      }else {
        Denominator <- Denominator + (1 - W) %*% exp(Ze_Normal[[k]] %*% theta)
      }
    }
    Denominator <- Denominator / r
    # 注意：后半部分Z1bar【不能】直接用colMeans！m1要除的是n1，不是N！
    HYJ_val <- Numerator / drop(Denominator) - Ze1_bar_al_Normal
    return(HYJ_val)
  }

}

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
  # ASMD.Denominator1 <- sqrt(colSums((t(t(Z1) - colMeans(Z1)))^2) / (nrow(Z1) - 1))

  ##ASMD
  ASMD <- abs(ASMD.Numerator_trt - ASMD.Numerator_col) / ASMD.Denominator1

  return(data.frame(X1ASMD = ASMD[1], X2ASMD = ASMD[2],
                    U3ASMD = ASMD[3], U4ASMD = ASMD[4]))
}

# ================================== #

{
  # Package Loading----
  library(MASS)
  library(WeightIt)
  library(reshape2)
  library(dplyr)
  library(gtools)
  library(ggplot2)
  library(stats)
  library(nleqslv)
  library(ebal)
  library(rootSolve)
  library(ggpubr)

  # Variable Information----
{
  ## Sample size
  N <- 10000
  # N <- 1000
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
  # N_Sim <- 500
  # N_Sim <- 10
}

  # Construct data frame for saving results----
{
  ## Normal distributed measurement errors
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
  # mat_SingleSim_oracle_Normal <- mat_SingleSim_eb_Normal
  ### mat_avgSim: [每一行]记录var_eps × rhoZ的某一个组合下，500次simulation结果取平均值后的结果
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
  # mat_avgSim_oracle_Normal <- mat_avgSim_eb_Normal
  ### mat_BiasSdMSE: [每一行]记录var_eps × rhoZ的某一个组合下，500次simulation后估计量的bias, sd, MSE.
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
  iterbfgsind.r2_Normal <- rep(TRUE, N_Sim)
  # iterbfgsind.r_Normal_rho06_sdeps06 <- iterbfgsind.r_Normal
  iterbfgsnum.r_Normal <- iterbfgsind.r_Normal
  iterbfgsnum.r2_Normal <- iterbfgsind.r_Normal
  ### BC - approximate bias corrector
  iterbcind.r_Normal <- iterbfgsind.r_Normal
  ### HCC
  iterhccind.r_Normal <- iterbfgsind.r_Normal
  iterhccpres.r_Normal <- iterbfgsind.r_Normal
  iterhccnum.r_Normal <- iterbfgsind.r_Normal
  ### HYJ
  iterhyjind.r_Normal <- iterbfgsind.r_Normal
  iterhyjpres.r_Normal <- iterbfgsind.r_Normal
  iterhyjnum.r_Normal <- iterbfgsind.r_Normal
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
  Count_Row <- 0
  ## The ith combination of 1-RR and rhoZ
  ith_vareps_rhoZ <- 0
}

# Correlation coefficient between X and U
rhoZ <- 0.3

# Add correlation to the covariance matrix of covariates
{
  Sigma_Z[1, 3] <- rhoZ
  Sigma_Z[3, 1] <- rhoZ
  Sigma_Z[2, 3] <- rhoZ
  Sigma_Z[3, 2] <- rhoZ
}

set.seed(seed = 1023)

for (var_eps in vec_Vareps) {
  # # Varaince of measurement error
  # var_eps <- 0.05

  # Covariance matrix of measurement errors containing only error-prone covariates.
  Sigma_e <- diag(rep(var_eps, 2))
  ## Sigma: the covariance matrix of measurement errorss containing both error-prone and error-free covariates.
  Sigma <- diag(c(diag(Sigma_e), 0, 0))

  # mat_SingleSim_ASMD_CEB: markdown ASMD in each simulation
  ## ncol = dim(Z)
  mat_SingleSim_ASMD_CEB <- matrix(0, nrow = N_Sim, ncol = 4)
  mat_SingleSim_ASMD_naiveEB <- matrix(0, nrow = N_Sim, ncol = 4)
  mat_SingleSim_ASMD_BCEB <- matrix(0, nrow = N_Sim, ncol = 4)
  mat_SingleSim_ASMD_HL <- matrix(0, nrow = N_Sim, ncol = 4)
  mat_SingleSim_ASMD_HW <- matrix(0, nrow = N_Sim, ncol = 4)

  for (i in 1:N_Sim) {
    # i <- 1

    # Mark the rhoZ
    ## Normal distributed measurement errors
    mat_SingleSim_ebme_Normal[i, "rhoZ"] <- rhoZ

    # Mark the Var_eps
    ## Normal distributed measurement errors
    mat_SingleSim_ebme_Normal[i, "Var_eps"] <- var_eps

    # True covariates
  {
    Z <- mvrnorm(n = N, mu = c(Mean_X, Mean_U), Sigma = Sigma_Z)
    trueX <- Z[, 1:p_X]
    U <- Z[, (p_X + 1):p_Z]
  }

    # Construct replicated covariate dataset
  {
    # Normally distributed measurement errors for replicated dataset
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

    # Estimate the variance-covariance matrix $\Sigma$
  {
    ## Normal distributed measurement errors
    meanZe_al_Normal <- Reduce("+", Ze_Normal) / r
    erhat_Normal <- matrix(rep(0, ncol(Z)^2), ncol = ncol(Z), nrow = ncol(Z))
    for (kk in 1:r) {
      erhat_Normal <- erhat_Normal +
        t(Ze_Normal[[kk]] - meanZe_al_Normal) %*% (Ze_Normal[[kk]] - meanZe_al_Normal)
    }
    Sigmahat_al_Normal <- erhat_Normal / ((r - 1) * N)
  }
    # round(Sigmahat_al_Normal, 3)

    # ZQY-EB: native Entropy Balancing
  {
    ## Normal distributed measurement errors
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
    # Replicate 1
    ## Normal distributed measurement errors
    opt.out.bfgs.r_Normal <- optim(par = rep(0, (ncol(Z))),
                                   fn = objective.ebme.r_Normal,
                                   gr = gradient.ebme.r_Normal,
                                   method = "BFGS",
                                   control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200))
    iterbfgsnum.r_Normal[i] <- opt.out.bfgs.r_Normal$counts[1]
    iterbfgsind.r_Normal[i] <- (max(gradient.ebme.r_Normal(opt.out.bfgs.r_Normal$par)) < 1e-04)
    ## double check
    if (is.na(iterbfgsind.r_Normal[i])) {
      iterbfgsind.r_Normal[i] <- FALSE
      print("Theta of CEB - NA clear!")
    }

    # -------------------- #

    # Replicate 2
    ## Normal distributed measurement errors
    opt.out.bfgs.r2_Normal <- optim(par = rep(0, (ncol(Z))),
                                    fn = objective.ebme.r2_Normal,
                                    gr = gradient.ebme.r2_Normal,
                                    method = "BFGS",
                                    control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200))
    iterbfgsnum.r2_Normal[i] <- opt.out.bfgs.r2_Normal$counts[1]
    iterbfgsind.r2_Normal[i] <- (max(gradient.ebme.r_Normal(opt.out.bfgs.r2_Normal$par)) < 1e-04)
    ## double check
    if (is.na(iterbfgsind.r2_Normal[i])) {
      iterbfgsind.r2_Normal[i] <- FALSE
      print("Theta of CEB - NA clear!")
    }

  }

    # BCEB
  {
    ## Normal distributed measurement errors
    # theta_bc_Normal <- drop(solve(hessian.zqy.r_Normal(theta0_Normal) - Sigmahat_al_Normal) %*%
    #                           hessian.zqy.r_Normal(theta0_Normal) %*%
    #                           theta0_Normal)
    theta_bc_Normal <- drop(solve(hessian.zqy.r_Normal(theta0_Normal) - Sigma) %*%
                              hessian.zqy.r_Normal(theta0_Normal) %*%
                              theta0_Normal)
    ### 判断BC是否无法求解：（默认【收敛】）存在以下情况任意其一都被判为【无解/不收敛】
    iterbcind.r_Normal[i] <- !(NaN %in% theta_bc_Normal | (sum(is.infinite(theta_bc_Normal)) > 0))
  }

    # HCC(2012)
  {
    ## Normal distributed measurement errors
    opt.out.hcc.r_Normal <- multiroot(f = gradient.hcc.r_Normal, start = rep(0, p_Z),
                                      rtol = .Machine$double.eps, maxiter = 200)
    iterhccind.r_Normal[i] <- (max(gradient.hcc.r_Normal(opt.out.hcc.r_Normal$root)) < sqrt(.Machine$double.eps)) & (opt.out.hcc.r_Normal$estim.precis != "NaN")
    iterhccpres.r_Normal[i] <- opt.out.hcc.r_Normal$estim.precis
    iterhccnum.r_Normal[i] <- opt.out.hcc.r_Normal$iter
  }

    # HYJ(2012)
  {
    ## Normal distributed measurement errors
    opt.out.hyj.r_Normal <- multiroot(f = gradient.hyj.r_Normal, start = rep(0, p_Z),
                                      rtol = .Machine$double.eps, maxiter = 200)
    iterhyjind.r_Normal[i] <- (max(gradient.hyj.r_Normal(opt.out.hyj.r_Normal$root)) < sqrt(.Machine$double.eps)) & (opt.out.hyj.r_Normal$estim.precis != "NaN")
    iterhyjpres.r_Normal[i] <- opt.out.hyj.r_Normal$estim.precis
    iterhyjnum.r_Normal[i] <- opt.out.hyj.r_Normal$iter
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
    ## Normal distributed measurement errors
    par_zqy_Normal <- as.matrix(opt.out.zqy.r_Normal$par, ncol = 1)
    weight_zqy_Normal <- exp(Ze0_Normal[[1]] %*% par_zqy_Normal)
    weight_zqy_Normal <- weight_zqy_Normal / sum(weight_zqy_Normal)
    ATT_zqy_Normal <- mean(Y1) - t(weight_zqy_Normal) %*% Y0
    mat_SingleSim_eb_Normal[i, "ATT"] <- mat_SingleSim_eb_Normal[i, "ATT"] + ATT_zqy_Normal
  }

    # CEB / EBME - bfgs C2: the proposed method
  {
    ## Normal distributed measurement errors
    # par_bfgs_Normal <- as.matrix(opt.out.bfgs.r_Normal$par, ncol = 1)
    par_bfgs_Normal <- list()
    par_bfgs_Normal[[1]] <- as.matrix(opt.out.bfgs.r_Normal$par, ncol = 1)
    par_bfgs_Normal[[2]] <- as.matrix(opt.out.bfgs.r2_Normal$par, ncol = 1)
    # weight_bfgs_Normal <- exp(Ze0_Normal[[1]] %*% par_bfgs_Normal)
    weight_bfgs_Normal <- matrix(rep(0, nrow(weight_zqy_Normal)))
    for (k in 1:r) {
      weight_bfgs_Normal <- weight_bfgs_Normal + exp(Ze0_Normal[[k]] %*% par_bfgs_Normal[[k]])
    }
    weight_bfgs_Normal <- weight_bfgs_Normal / sum(weight_bfgs_Normal)
    ATT_bfgs_Normal <- mean(Y1) - t(weight_bfgs_Normal) %*% Y0
    mat_SingleSim_ebme_Normal[i, "ATT"] <- mat_SingleSim_ebme_Normal[i, "ATT"] + ATT_bfgs_Normal
  }

    # BCEB
  {
    ## Normal distributed measurement errors
    par_bc_Normal <- as.matrix(theta_bc_Normal, ncol = 1)
    weight_bc_Normal <- exp(Ze0_Normal[[1]] %*% par_bc_Normal)
    weight_bc_Normal <- weight_bc_Normal / sum(weight_bc_Normal)
    ATT_bc_Normal <- mean(Y1) - t(weight_bc_Normal) %*% Y0
    mat_SingleSim_bc_Normal[i, "ATT"] <- mat_SingleSim_bc_Normal[i, "ATT"] + ATT_bc_Normal
  }

    # HCC
  {
    ## Normal distributed measurement errors
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
    ## Normal distributed measurement errors
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

    # mat_SingleSim_ASMD_CEB
    mat_SingleSim_ASMD_CEB[i,] <- round(unlist(ASMDCal(weight_bfgs_Normal)), 3)
    mat_SingleSim_ASMD_naiveEB[i,] <- round(unlist(ASMDCal(weight_zqy_Normal)), 3)
    mat_SingleSim_ASMD_BCEB[i,] <- round(unlist(ASMDCal(weight_bc_Normal)), 3)
    mat_SingleSim_ASMD_HL[i,] <- round(unlist(ASMDCal(weight_hcc_Normal)), 3)
    mat_SingleSim_ASMD_HW[i,] <- round(unlist(ASMDCal(weight_hyj_Normal)), 3)

    print(paste0(i, "/", N_Sim, ", the ", i, "th Simulation for var_eps = ", var_eps, " finished."))
  }

  # mat_ASMD_CEB_Normal
  mat_ASMD_CEB_Normal <- rbind(mat_ASMD_CEB_Normal,
                               colMeans(mat_SingleSim_ASMD_CEB, na.rm = TRUE))
  mat_ASMD_naiveEB_Normal <- rbind(mat_ASMD_naiveEB_Normal,
                                   colMeans(mat_SingleSim_ASMD_naiveEB, na.rm = TRUE))
  mat_ASMD_BCEB_Normal <- rbind(mat_ASMD_BCEB_Normal,
                                colMeans(mat_SingleSim_ASMD_BCEB, na.rm = TRUE))
  mat_ASMD_HL_Normal <- rbind(mat_ASMD_HL_Normal,
                              colMeans(mat_SingleSim_ASMD_HL, na.rm = TRUE))
  mat_ASMD_HW_Normal <- rbind(mat_ASMD_HW_Normal,
                              colMeans(mat_SingleSim_ASMD_HW, na.rm = TRUE))
  # round(mat_ASMD_CEB_Normal, 3)

  # 矩阵清空为下一轮simulation做准备
  mat_SingleSim_ASMD_CEB <- matrix(0, nrow = N_Sim, ncol = 4)

}

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

round(mat_ASMD_CEB_Normal, 3)
round(mat_ASMD_naiveEB_Normal, 3)
round(mat_ASMD_BCEB_Normal, 3)
round(mat_ASMD_HL_Normal, 3)
round(mat_ASMD_HW_Normal, 3)

# mat_ASMD_CEB_Normal
# mat_ASMD_naiveEB_Normal
# mat_ASMD_BCEB_Normal
# mat_ASMD_HL_Normal
# mat_ASMD_HW_Normal