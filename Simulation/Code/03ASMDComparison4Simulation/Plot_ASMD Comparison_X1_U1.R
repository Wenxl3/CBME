# -*- coding: utf-8 -*-
# ------------------------------------------
# Title       : ASMD Comparison_X1_U1
# Objective   : Draw the plot of ASMD Comparison.
# Remark - 1  : 2 fixed var_eps scenarios: (1)0.2; (2)0.5.
# Remark - 2  : Measurement error follows the normal distribution.
# Remark - 3  : x-lab: N \in {1k, 1.5k, 2k, ..., 1w}.
# Remark - 4  : y-lab: Compare 5 methods: naive EB, CEB, BCEB, HL, HW.
# Remark - 5  : Only compare the ASMD of X1 and U1.
# Remark - 6  : Four subplots.
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
  ##  ASMD Calculation Function----
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
  source(file = "./CEB/R/03ASMDComparison4Simulation/ObjectiveGradient_NormalME.R")

  # Variable Information----
{
  ## Sample size
  vec_N <- c(500, seq(1000, 10000, by = 1000))
  ## numbers of replicate, i.e. N_Replicate
  r <- 2
  ## Covariates----
  ### X is error-prone data
  Mean_X <- c(4, 2)
  Mean_U <- c(3, 1)
  p_X <- length(Mean_X)
  p_U <- length(Mean_U)
  p_Z <- p_U + p_X
  ### Sigma_Z: Covariance matrix of covariates
  Sigma_Z <- diag(rep(1, 4))
  ### rhoZ: Correlation coefficient between error-prone covariates and error-free covariates
  # vec_rhoZ <- 3 / 10
  ### Variance of measurement errors
  # vec_Vareps <- seq(0, 0.5, by = 0.05)
  vec_Vareps <- c(0.2, 0.5)
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

  # Scenarios Count
  TotalScenCount <- length(vec_N) * length(vec_Vareps)
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
  ## Normal distributed measurement errors
{
  ## Iteration numbers and indicator pool
  ### EBME-BFGS
  itercebind.r_Normal <- rep(TRUE, N_Sim)
  ### BC - approximate bias corrector
  iterbcind.r_Normal <- itercebind.r_Normal
  ### HCC
  iterhccind.r_Normal <- itercebind.r_Normal
  ### HYJ
  iterhyjind.r_Normal <- itercebind.r_Normal

  # matrix markdown the ASMD
  # 6 cols: var_eps, N, avgASMD_X1, avgASMD_X2, avgASMD_U1, avgASMD_U2
  mat_ASMD_naiveEB_Normal <- NULL
  mat_ASMD_CEB_Normal <- NULL
  mat_ASMD_BCEB_Normal <- NULL
  mat_ASMD_HL_Normal <- NULL
  mat_ASMD_HW_Normal <- NULL
}

}

  # Simulation progress counting----
  pg <- 0
  ## The ith combination of 1-RR and rhoZ
  ith_N_vareps_rhoZ <- 0
  mat_ASMD <- NULL
}

# Simulation----
{
  rhoZ <- 0.3

  # Add correlation to the covariance matrix of covariates
  Sigma_Z[1, 3] <- rhoZ
  Sigma_Z[3, 1] <- rhoZ
  Sigma_Z[2, 3] <- rhoZ
  Sigma_Z[3, 2] <- rhoZ

  # mat_SingleSim_ASMD_CEB: markdown ASMD in each simulation
  ## ncol = dim(Z)
  mat_SingleSim_ASMD_naiveEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
  mat_SingleSim_ASMD_CEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
  mat_SingleSim_ASMD_BCEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
  mat_SingleSim_ASMD_HL <- matrix(0, nrow = N_Sim, ncol = p_Z)
  mat_SingleSim_ASMD_HW <- matrix(0, nrow = N_Sim, ncol = p_Z)
}

for (var_eps in vec_Vareps) {
  # # Varaince of measurement error
  # var_eps <- 0.5

  # Covariance matrix of measurement errors containing only error-prone covariates.
  Sigma_e <- diag(rep(var_eps, p_X))
  ## Sigma: the covariance matrix of measurement errorss containing both error-prone and error-free covariates.
  Sigma <- diag(c(diag(Sigma_e), rep(0, p_U)))

  for (N in vec_N) {
    # N <- 1000

    for (i in 1:N_Sim) {
      # i <- 1

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
      ## Normal distributed measurement errors
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
      # theta_CEB - solve the objective funtion
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
        itercebind.r_Normal[i] <- (max(abs(gradient.ebme.r_Normal(theta_CEB_solveObj))) < 1e-04)
        ## double check
        if (is.na(itercebind.r_Normal[i]) | is.nan(itercebind.r_Normal[i])) {
          itercebind.r_Normal[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }else {
        theta_CEB_Normal <- theta_CEB_solveGradient
        itercebind.r_Normal[i] <- (max(abs(output.CEB.solveGradient_Normal$f.root)) < sqrt(.Machine$double.eps)) & (output.CEB.solveGradient_Normal$estim.precis != "NaN")
        ## double check
        if (is.na(itercebind.r_Normal[i]) | is.nan(itercebind.r_Normal[i])) {
          itercebind.r_Normal[i] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }
    }

      # BCEB
    {
      ## Normal distributed measurement errors
      theta_bc_Normal <- drop(solve(hessian.zqy.r_Normal(theta0_Normal) - Sigma) %*%
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
      ## Normal distributed measurement errors
      weight_zqy_Normal <- exp(Ze0_Normal[[1]] %*% theta0_Normal)
      weight_zqy_Normal <- weight_zqy_Normal / sum(weight_zqy_Normal)
    }

      # CEB / EBME - bfgs C2: the proposed method
    {
      ## Normal distributed measurement errors
      weight_bfgs_Normal <- exp(Ze0_Normal[[1]] %*% theta_CEB_Normal)
      weight_bfgs_Normal <- weight_bfgs_Normal / sum(weight_bfgs_Normal)
    }

      # BCEB
    {
      ## Normal distributed measurement errors
      weight_bc_Normal <- exp(Ze0_Normal[[1]] %*% theta_bc_Normal)
      weight_bc_Normal <- weight_bc_Normal / sum(weight_bc_Normal)
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
    }
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
      print(paste0("rhoZ = ", rhoZ, ", var_eps = ", var_eps, ", N = ", N, ", Sim.i = ", i))
      print(paste0("Progress = ",
                   round(pg / (length(vec_Vareps) *
                     length(vec_N) *
                     N_Sim) * 100,
                         digits = 2), "%"))
    }

    ith_N_vareps_rhoZ <- ith_N_vareps_rhoZ + 1
    print(paste0("The ", ith_N_vareps_rhoZ, "-th(/", TotalScenCount, ") row of mat_avgSim value assignment -- start"))

    # matrix of average value of ASMD
  {
    mat_ASMD_naiveEB_Normal <- rbind(mat_ASMD_naiveEB_Normal,
                                     c(var_eps, N, colMeans(mat_SingleSim_ASMD_naiveEB, na.rm = TRUE)))
    mat_ASMD_CEB_Normal <- rbind(mat_ASMD_CEB_Normal,
                                 c(var_eps, N, colMeans(mat_SingleSim_ASMD_CEB[itercebind.r_Normal,], na.rm = TRUE)))
    mat_ASMD_BCEB_Normal <- rbind(mat_ASMD_BCEB_Normal,
                                  c(var_eps, N, colMeans(mat_SingleSim_ASMD_BCEB[iterbcind.r_Normal,], na.rm = TRUE)))
    mat_ASMD_HL_Normal <- rbind(mat_ASMD_HL_Normal,
                                c(var_eps, N, colMeans(mat_SingleSim_ASMD_HL[iterhccind.r_Normal,], na.rm = TRUE)))
    mat_ASMD_HW_Normal <- rbind(mat_ASMD_HW_Normal,
                                c(var_eps, N, colMeans(mat_SingleSim_ASMD_HW[iterhyjind.r_Normal,], na.rm = TRUE)))
    # round(mat_ASMD_CEB_Normal, 3)
  }


  {
    ## EBME-BFGS
    itercebind.r_Normal <- rep(TRUE, N_Sim)
    ## BC - approximate bias corrector
    iterbcind.r_Normal <- itercebind.r_Normal
    ## HCC
    iterhccind.r_Normal <- itercebind.r_Normal
    ## HYJ
    iterhyjind.r_Normal <- itercebind.r_Normal

    #
    mat_SingleSim_ASMD_CEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_naiveEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_BCEB <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HL <- matrix(0, nrow = N_Sim, ncol = p_Z)
    mat_SingleSim_ASMD_HW <- matrix(0, nrow = N_Sim, ncol = p_Z)

  }

    # end  for (N in vec_N)
  }

  # end for (var_eps in vec_Vareps)
}

# Save .Rdata
save.image(file = "Plot_ASMDComparison_NormalEps_Original.RData")


# Plot output----
{
  ## mat_ASMD
  mat_ASMD_eb_Normal <- cbind(mat_ASMD_naiveEB_Normal, Method = "Naive")
  mat_ASMD_ebme_Normal <- cbind(mat_ASMD_CEB_Normal, Method = "CEB")
  mat_ASMD_bc_Normal <- cbind(mat_ASMD_BCEB_Normal, Method = "BCEB")
  mat_ASMD_hcc_Normal <- cbind(mat_ASMD_HL_Normal, Method = "HL")
  mat_ASMD_hyj_Normal <- cbind(mat_ASMD_HW_Normal, Method = "HW")
  avgmat_ASMD_Normal <- rbind(mat_ASMD_eb_Normal,
                              mat_ASMD_ebme_Normal,
                              mat_ASMD_bc_Normal,
                              mat_ASMD_hcc_Normal,
                              mat_ASMD_hyj_Normal)
  colnames(avgmat_ASMD_Normal) <- c("var_eps", "N", paste0("X", 1:p_X), paste0("U", 1:p_U), "Method")
  csvName_mat_ASMD <- paste0("mat_NormalEps_avgASMD_rhoZ", rhoZ, "_Sim", N_Sim,
                             ".xlsx")
  openxlsx::write.xlsx(avgmat_ASMD_Normal,
                       file = csvName_mat_ASMD,
                       row.names = FALSE)
}

# Plots----
{
  # Only save the ASMD data of X1 and U1
  # library(dplyr)
  ## 'select' is only available for data.frame variable
  colnames_rm <- paste0(c("X", "U"), 2)
  avgmat_ASMD_Normal_X1_U1 <- select(as.data.frame(avgmat_ASMD_Normal), -colnames_rm)

  ASMDdf_Normal <- melt(as.data.frame(avgmat_ASMD_Normal_X1_U1),
                        id.vars = c("var_eps", "N", "Method"))
  ASMDdf_Normal$value <- as.numeric(ASMDdf_Normal$value)
  which(is.na(ASMDdf_Normal$value))
  ASMDdf_Normal$Method <- factor(ASMDdf_Normal$Method, levels = c("Naive", "CEB", "BCEB", "HL", "HW"))
  ASMDdf_Normal$N <- as.numeric(ASMDdf_Normal$N)

  rownindex_vare02 <- which(ASMDdf_Normal$var_eps == 0.2)
  rownindex_vare05 <- which(ASMDdf_Normal$var_eps == 0.5)
  ASMDdf_Normal[rownindex_vare02, "var_eps"] <- "Var = 0.2"
  ASMDdf_Normal[rownindex_vare05, "var_eps"] <- "Var = 0.5"
  ASMDdf_Normal$var_eps <- as.factor(ASMDdf_Normal$var_eps)

  p_ASMDdf_X1_U1_NormalEps <- ggplot(ASMDdf_Normal, aes(x = N, y = value, group = Method)) +
    facet_grid(rows = vars(variable), cols = vars(var_eps),
               scales = "free_y") +
    geom_line(aes(x = N, y = value, color = Method, group = Method)) +
    geom_point(aes(x = N, y = value, color = Method, shape = Method)) +
    labs(x = "Sample Size", y = "Absolute Standardized Mean Difference (ASMD)", title = "") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    scale_x_continuous(breaks = c(500, seq(1000, 10000, by = 1000)), labels = c(500, seq(1000, 10000, by = 1000))) +
    # ylim(-0.8, 0.6) +
    scale_y_continuous(breaks = seq(0, max(ASMDdf_Normal$value) + 0.03, by = 0.01), labels = seq(0, max(ASMDdf_Normal$value) + 0.03, by = 0.01)) +
    theme_bw()
  p_ASMDdf_X1_U1_NormalEps
  ggsave(filename = "p_ASMD_X1_U1_NormalEps.pdf", p_ASMDdf_X1_U1_NormalEps,
         width = 26, height = 23,
         units = "cm")
}

# Plots without Naive EB----
{
  ASMDdf_Normal_NoNaiveEB <- ASMDdf_Normal
  rowindex_naiveEB <- which(ASMDdf_Normal_NoNaiveEB$Method == "Naive")
  ASMDdf_Normal_NoNaiveEB <- ASMDdf_Normal_NoNaiveEB[-rowindex_naiveEB,]

  p_ASMDdf_X1_U1_NormalEps_NoNaiveEB <- ggplot(ASMDdf_Normal_NoNaiveEB, aes(x = N, y = value, group = Method)) +
    facet_grid(rows = vars(variable), cols = vars(var_eps),
               scales = "free_y") +
    geom_line(aes(x = N, y = value, color = Method, group = Method)) +
    geom_point(aes(x = N, y = value, color = Method, shape = Method)) +
    labs(x = "Sample Size", y = "Absolute Standardized Mean Difference (ASMD)", title = "") +
    geom_abline(slope = 0, intercept = 0, lty = 3, lwd = 0.5) +
    scale_x_continuous(breaks = seq(0, 11000, 1000)) +
    theme_bw()
  p_ASMDdf_X1_U1_NormalEps_NoNaiveEB
  ggsave(filename = "p_ASMDdf_X1_U1_NormalEps_NoNaiveEB.pdf", p_ASMDdf_X1_U1_NormalEps_NoNaiveEB,
         width = 26, height = 23,
         units = "cm")
}

# Save .Rdata
save.image(file = "Compare_avgASMD_X1_U1_NormalEps.RData")

