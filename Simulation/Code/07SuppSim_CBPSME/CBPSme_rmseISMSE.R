# -*- coding: utf-8 -*-
# ------------------------------------------
# Title       : 17.2.KSdata-4Covs_ATE_Y-a_PS-a
# Objective   : Misspecification model + GMM solution
## Misspecification model: Y_scen = "a", PS_scen = "a"
# Code Reference(s): A framework for covariate balance using Bregman distances
# Created on  : 2023/03/22
# ------------------------------------------

{
  # ====== Basic Setting ======= #
{
  # Working Directory
  rm(list = ls())
  getwd()
  # "warn" off
  options(warn = -1)

  Time.startpoint <- Sys.time()

  # Package
  ## Multivariate Normal Distribution
  library(MASS)
  ## CBPS
  library(CBPS)
  # ## Solve the nonlinear equation (for example, for DiShu)
  # library(nleqslv)
  ## Create the progress bar
  library(progress)
  ## reshaple the data frame
  library(reshape2)
  ## parallel
  library(parallel)
  library(pbapply)
  ## compiler: accelerate function calculation
  library(compiler)
}


  # ====== Self-designed Function ======= #

{
  # =========================================================== #
  # ****** -------------------- RMSE ------------------- ****** #
  # =========================================================== #

  RMSEfun <- function(actual, predicted,
                      NA.rm = T) {
    SE <- (actual - predicted)^2
    MSE <- mean(SE, na.rm = NA.rm)
    # RMSE <- sqrt(MSE)

    return(MSE)
  }

  # =========================================================== #
  # ****** --- UniSimDataset construction 4 parallel --- ****** #
  # =========================================================== #

  UniSimDataset <- function(ith_sim,
                            N = N,
                            var_eps = varE) {

    # var_eps <- 0.1

    # load R packages
    library(MASS)

    # ****** ------ Kang & Schafer (2007) simulate scenarios ------ ****** #

    # simulate scenarios
    ks_data <- function(tau, n, sig2, rho, y_scen = c("a", "b"), z_scen = c("a", "b")) {

      # covariates
      x1 <- stats::rnorm(n, 0, 1)
      x2 <- stats::rnorm(n, 0, 1)
      x3 <- stats::rnorm(n, 0, 1)
      x4 <- stats::rnorm(n, 0, 1)

      # transformed predictors
      u1 <- as.numeric(scale(exp(x1 / 2)))
      u2 <- as.numeric(scale(x2 / (1 + exp(x1)) + 10))
      u3 <- as.numeric(scale((x1 * x3 / 25 + 0.6)^3))
      u4 <- as.numeric(scale((x2 + x4 + 20)^2))

      # treatment probabilities
      if (z_scen == "b")
        e_X <- 1 / (1 + exp(-(-u1 + 0.5 * u2 - 0.25 * u3 - 0.1 * u4)))
      else
        e_X <- 1 / (1 + exp(-(-x1 + 0.5 * x2 - 0.25 * x3 - 0.1 * x4)))

      r_exp <- stats::runif(n)
      z <- ifelse(r_exp < e_X, 1, 0)

      # error variance
      R <- matrix(rho, nrow = 2, ncol = 2)
      diag(R) <- 1
      V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
      Sig <- V %*% R %*% V

      if (y_scen == "b")
        mu <- 210 +
          27.4 * u1 +
          13.7 * u2 +
          13.7 * u3 +
          13.7 * u4
      else
        mu <- 210 +
          27.4 * x1 +
          13.7 * x2 +
          13.7 * x3 +
          13.7 * x4

      eval <- eigen(Sig, symmetric = TRUE)
      y_init <- matrix(stats::rnorm(n * 2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
      y_tmp <- t(eval$vectors %*%
                   diag(sqrt(eval$values), nrow = 2) %*%
                   t(y_init)) # SVD
      y_pot <- y_tmp + cbind(mu, mu + tau) # include causal effect

      # observed outcome
      y <- z * y_pot[, 2] + (1 - z) * y_pot[, 1]

      # create simulation dataset
      sim.dataset <- as.data.frame(cbind(y, z, x1, x2, x3, x4, u1, u2, u3, u4))

      sim.data <- list(ks_data = sim.dataset,
                       ps = e_X)

      return(sim.data)

    }

    # ====== Global Parameter Setting ====== #
  {
    dim_X <- 2
    dim_V <- 2
    dim_Z <- dim_X + dim_V

    # Outcome variable
    ## Covariance matrix of potential outcomes
    Sigma_Y <- 10
    ## Correlation coefficient matrix of potential outcomes
    rhoY01 <- 0
    ## True value of average treatment effect
    TrueValue_ATE <- 10
    ## True value of average effect of population
    TrueValue_mu <- 210


    # Set the scenario of outcome regression model & propensity score model
    ## "a": use {x1, x2} as predictors
    ## "b": use {u1, u2} as predictors
    Y_scen <- "a"
    PS_scen <- "a"
  }

    set.seed(seed = ith_sim)

    UniSimData <- ks_data(tau = TrueValue_ATE,
                          n = N,
                          sig2 = Sigma_Y, rho = rhoY01,
                          y_scen = Y_scen, z_scen = PS_scen)
    Dataset <- UniSimData$ks_data

    Predictors <- Dataset
    Predictors$y <- NULL
    Predictors$z <- NULL

    # Observed dataset without measurement error
    ## We can only obsever the "x" covariates no matter the scenario of PS or Y
    Dataset_obs_withoutME <- Dataset[, c("y", "z", "x1", "x2", "x3", "x4")]
    ## X1, X2: error-porn covariates
    ## V1, V2: error-free covariates
    colnames(Dataset_obs_withoutME) <- c("Y", "W", "X1", "X2", "V1", "V2")
    X1 <- Dataset_obs_withoutME$X1
    X2 <- Dataset_obs_withoutME$X2
    V1 <- Dataset_obs_withoutME$V1
    V2 <- Dataset_obs_withoutME$V2
    Covs_X <- cbind(X1, X2)
    Covs_V <- cbind(V1, V2)
    # Z <- cbind(X1, X2, V1, V2)
    Z <- cbind(Covs_X, Covs_V)
    Y <- Dataset_obs_withoutME$Y
    # W: the treatment indicator
    W <- Dataset_obs_withoutME$W
    TrtIndex <- which(W == 1)
    ColIndex <- which(W == 0)
    n1 <- length(TrtIndex)
    n0 <- length(ColIndex)

    # Measurement error (ME)
    ifelse(dim_X == 1,
           Sigma_e <- var_eps,
           Sigma_e <- diag(x = rep(var_eps, dim_X)))
    ME <- mvrnorm(n = N,
                  mu = rep(0, dim_X),
                  Sigma = Sigma_e)
    # ME <- matrix(data = ME, nrow = N, ncol = dim_X)
    Covs_Xstar <- Covs_X + ME
    Zstar <- data.frame(cbind(Covs_Xstar, Covs_V))
    # colnames(Zstar) <- c("X1", "X2", "V1", "V2")

    # Treatment and Propensity Score
    ps <- UniSimData$ps

    # Observed dataset with measurement error
    Dataset_obs <- cbind(Y, W, Zstar)
    colnames(Dataset_obs) <- c("Y", "W", paste0("X", 1:2), paste0("V", 1:2))

    Data_obs <- list(Predictors = Predictors,
                     Dataset_obs = Dataset_obs,
                     ps = ps)

    return(Data_obs)
  }


  # ================================================ #
  # ****** ------ UniSimATE 4 parallel ------ ****** #
  # ================================================ #

  UniSimATE <- function(data,
                        N = N,
                        var_eps = varE) {

    # var_eps <- 0.1
    # data <- UniSimDataset(ith_sim = 1, var_eps = 0.1)
    # data <- AllSimData[[1]]

    # ====== load R packages ====== #

    library(CBPS)
    library(compiler)

    # ====== Self-write Function ====== #

  {
    # ****** ------ optim with compiler ------ ****** #
    optim.compiler <- cmpfun(f = optim)
    optimize.compiler <- cmpfun(f = optimize)

    # ****** ------ JICBPS ------ ****** #
    ##Loss function for balance constraints, returns the squared imbalance along each dimension.

    GMM.JICBPS <- function(para_beta,
                           x) {
      # para_beta <- c(para_1 = 0.2, para_X = 0.5, para_V = 1)
      # x <- unisimdataset

      # There is probably an X'X missing somewhere that is causing these variance problems.
      probs.min <- 1e-6

      Covs <- unisimdataset
      Covs$Y <- NULL
      Covs$W <- NULL
      Covs <- as.matrix(cbind(1, Covs))

      para_beta <- matrix(para_beta)

      theta.curr <- as.vector(Covs %*% as.matrix(para_beta))

      probs.curr <- (1 + exp(-theta.curr))^-1
      probs.curr <- pmin(1 - probs.min, probs.curr)
      probs.curr <- pmax(probs.min, probs.curr)
      probs.curr <- as.vector(probs.curr)

      pi_beta <- probs.curr

      weights_JICBPS <- (W - pi_beta) / pi_beta / (1 - pi_beta)

      # A-function
      JICBPS_val <- Covs * weights_JICBPS

      return(JICBPS_val)
    }

    # ****** ------ OICBPS ------ ****** #
    # use "CBPS" original code

    GMM.OICBPS <- function(para_beta,
                           x) {
      # para_beta <- c(para_1 = 0.2, para_X = 0.5, para_V = 1)
      # x <- unisimdataset

      # There is probably an X'X missing somewhere that is causing these variance problems.
      probs.min <- 1e-6

      Covs <- unisimdataset
      Covs$Y <- NULL
      Covs$W <- NULL
      Covs <- as.matrix(cbind(1, Covs))

      para_beta <- matrix(para_beta)

      theta.curr <- as.vector(Covs %*% as.matrix(para_beta))

      probs.curr <- (1 + exp(-theta.curr))^-1
      probs.curr <- pmin(1 - probs.min, probs.curr)
      probs.curr <- pmax(probs.min, probs.curr)
      probs.curr <- as.vector(probs.curr)

      # weights 4 estimate ATE
      w.curr <- (probs.curr - 1 + W)^-1

      # Generate the of mean imbalance by weights for each covariate. - JICBPS value
      w.curr.del <- Covs * w.curr
      # GMM.JICBPS(para_beta = para_beta, x = x)

      # Generate score function value
      score.val <- Covs * (W - probs.curr)
      # sum(score.val[, 1] == Covs[, 2] * (W - probs.curr))
      # sum(score.val[, 2] == Covs[, 3] * (W - probs.curr))

      # Generate g-bar, as in the paper "CBPS".
      gbar <- cbind(score.val, w.curr.del)

      return(gbar)
    }

    # ****** ------ Hstar ------ ****** #

    Hstar <- function(para_beta, parms) {
      # para_beta <- c(para_1 = 0,
      #                para_X1 = -1, para_X2 = 0.5,
      #                para_V1 = -0.25, para_V2 = -0.1)
      # parms <- list(Xstar = Xstar, V = V, sigmaE = sigmaE)
      # parms <- list(Xstar = Xstar, V = V, var_eps = var_eps)

      # para_beta
      para_X1 <- para_beta["para_X1"]
      para_X2 <- para_beta["para_X2"]
      para_X <- rbind(para_X1, para_X2)

      Xstar <- as.matrix(parms[["Xstar"]])
      V <- parms[["V"]]
      # sigmaE <- parms[["sigmaE"]]
      # var_eps <- sigmaE^2

      # The measurement error covariance matrix containing only error-prone covariate.
      ifelse(dim_X == 1,
             Sigma_e <- var_eps,
             Sigma_e <- diag(x = rep(var_eps, dim_X)))

      delta_H <- Xstar - matrix(data = as.vector(Sigma_e %*% para_X * 0.5),
                                byrow = T,
                                nrow = nrow(Xstar), ncol = ncol(Xstar))

      CovInfo_H <- as.matrix(cbind(One = 1, delta_H = delta_H, V = V))
      Hstar_val <- exp(CovInfo_H %*% para_beta)

      return(Hstar_val)
    }

    # ****** ------ Sstar 4 DiShu solution ------ ****** #

    # S_star function for 【individual】 used for GMM
    GMM.Sstar <- function(para_beta,
                          x) {
      # para_beta <- c(para_1 = 0,
      #                para_X1 = -1, para_X2 = 0.5,
      #                para_V1 = -0.25, para_V2 = -0.1)
      # x <- unisimdataset

      # para_beta
      para_X1 <- para_beta["para_X1"]
      para_X2 <- para_beta["para_X2"]
      para_X <- rbind(para_X1, para_X2)
      para_beta <- as.matrix(para_beta)
      # dataset
      dataset <- as.data.frame(x)
      Xstar <- dataset[, c("X1", "X2")]
      V <- dataset[, c("V1", "V2")]
      # The measurement error covariance matrix containing only error-prone covariate.
      ifelse(dim_X == 1,
             Sigma_e <- var_eps,
             Sigma_e <- diag(x = rep(var_eps, dim_X)))

      # S-function
      delta <- Xstar + matrix(data = rep(Sigma_e %*% para_X, each = N), nrow = N) * (W - 0.5)
      CovInfo <- as.matrix(cbind(One = 1, delta = delta, V = V))
      temp <- W - 1 / (1 + exp(CovInfo %*% (-para_beta)))
      temp <- as.vector(temp)
      # colSums(cbind(temp, temp * V, temp * delta))
      # S_star <- abs(as.vector(t(temp) %*% CovInfo))
      # S_star <- cbind(temp, temp, temp, temp, temp) * CovInfo
      S_star <- temp * CovInfo

      return(S_star)
    }

    # ****** ------ Huang 4 Huang (2001) solution ------ ****** #

    #  objective.Huang2001 function for 【individual】 used for GMM
    GMM.Huang <- function(para_beta,
                          x) {

      # para_beta <- c(para_1 = 0,
      #                para_X1 = -1, para_X2 = 0.5,
      #                para_V1 = -0.25, para_V2 = -0.1)
      # x <- unisimdataset

      # para_beta
      para_X1 <- para_beta["para_X1"]
      para_X2 <- para_beta["para_X2"]
      para_X <- rbind(para_X1, para_X2)
      # para_beta <- as.matrix(para_beta)

      # dataset
      dataset <- as.data.frame(x)
      Xstar <- dataset[, c("X1", "X2")]
      V <- dataset[, c("V1", "V2")]

      parms <- list(Xstar = Xstar,
                    V = V,
                    var_eps = var_eps)

      # The measurement error covariance matrix containing only error-prone covariate.
      ifelse(dim_X == 1,
             Sigma_e <- var_eps,
             Sigma_e <- diag(x = rep(var_eps, dim_X)))

      # CovInfo
      CovInfo <- as.matrix(cbind(One = 1,
                                 Xstar = Xstar,
                                 V = V))
      # ------------------------------------------ #
      CovInfo2_phiminus <- as.matrix(cbind(One = 1,
                                           Xstar = Xstar + matrix(data = rep(Sigma_e %*% para_X, N), byrow = T, ncol = ncol(Xstar)),
                                           V = V))
      temp_minus <- W *
        as.vector(Hstar(para_beta = -para_beta, parms = parms)) *
        CovInfo2_phiminus
      # mean_Phi_minus_value <- colMeans((W - 1) * CovInfo + temp_minus)
      Phi_minus_value <- (W - 1) * CovInfo + temp_minus

      # ------------------------------------------ #
      CovInfo2_phiplus <- as.matrix(cbind(One = 1,
                                          Xstar = Xstar - matrix(data = rep(Sigma_e %*% para_X, N), byrow = T, ncol = ncol(Xstar)),
                                          V = V))
      temp_plus <- (W - 1) *
        as.vector(Hstar(para_beta = para_beta, parms = parms)) *
        CovInfo2_phiplus
      # mean_Phi_plus_value <- colMeans(W * CovInfo + temp_plus)
      Phi_plus_value <- W * CovInfo + temp_plus

      # mean_PhiVector <- matrix(data = c(mean_Phi_minus_value, mean_Phi_plus_value),
      #                          ncol = 1)
      PhiVector <- cbind(Phi_minus_value, Phi_plus_value)

      return(PhiVector)
    }

    # ****** ------ Astar 4 just-identified cCBPS solution ------ ****** #

    # A_star function for 【individual】 used for GMM
    GMM.Astar <- function(para_beta,
                          x) {
      # para_beta <- c(para_1 = 0,
      #                para_X1 = -1, para_X2 = 0.5,
      #                para_V1 = -0.25, para_V2 = -0.1)
      # x <- unisimdataset

      # para_beta
      para_beta.orig <- para_beta
      # para_beta
      para_X1 <- para_beta["para_X1"]
      para_X2 <- para_beta["para_X2"]
      para_X <- rbind(para_X1, para_X2)
      # para_beta <- as.matrix(para_beta)

      # dataset
      dataset <- as.data.frame(x)
      Xstar <- dataset[, c("X1", "X2")]
      V <- dataset[, c("V1", "V2")]

      parms <- list(Xstar = Xstar,
                    V = V,
                    var_eps = var_eps)

      # The measurement error covariance matrix containing only error-prone covariate.
      ifelse(dim_X == 1,
             Sigma_e <- var_eps,
             Sigma_e <- diag(x = rep(var_eps, dim_X)))

      # A-function
      c1 <- (W * (1 + Hstar(para_beta = -para_beta.orig, parms = parms)) +
        (W - 1) * (1 + Hstar(para_beta = para_beta.orig, parms = parms)))
      part1 <- as.matrix(cbind(c1, c1 * Xstar, c1 * V))
      c2 <- W * Hstar(para_beta = -para_beta.orig, parms = parms) -
        (W - 1) * Hstar(para_beta = para_beta.orig, parms = parms)
      # part2 <- as.matrix(cbind(0, c2 * sigmaE^2 * para_X, 0))
      part2 <- as.matrix(cbind(0,
                               matrix(data = rep(Sigma_e %*% para_X, N), byrow = T, ncol = ncol(Xstar)) * as.vector(c2),
                               matrix(data = 0, nrow = nrow(V), ncol = ncol(V))))
      # JIcCBPS_val
      # A_star <- abs(as.matrix(colSums(part1 + part2)))
      A_star <- part1 + part2

      return(A_star)
    }

    # ****** ------ OIcCBPS 4 over-identified cCBPS solution ------ ****** #

    # Overidentified cCBPS function for 【individual】 used for GMM
    GMM.OIcCBPS <- function(para_beta0,
                            x) {
      # para_beta0 <- c(para_1 = 0,
      #                para_X1 = -1, para_X2 = 0.5,
      #                para_V1 = -0.25, para_V2 = -0.1)
      # x <- unisimdataset

      # para_beta
      para_X1 <- para_beta0["para_X1"]
      para_X2 <- para_beta0["para_X2"]
      para_X <- rbind(para_X1, para_X2)

      # dataset
      dataset <- as.data.frame(x)
      Xstar <- dataset[, c("X1", "X2")]
      V <- dataset[, c("V1", "V2")]

      # The measurement error covariance matrix containing only error-prone covariate.
      ifelse(dim_X == 1,
             Sigma_e <- var_eps,
             Sigma_e <- diag(x = rep(var_eps, dim_X)))

      # S-function
      S_star <- GMM.Sstar(para_beta = para_beta0,
                          x = dataset)

      # A-function
      A_star <- GMM.Astar(para_beta = para_beta0,
                          x = dataset)

      # OIcCBPS_val
      OIcCBPS_val <- cbind(S_star, A_star)

      return(OIcCBPS_val)
    }

    # ****** ------ propensity score ------ ****** #

    # The propensity score estimating function
    # Input: the estimating coefficients of propensity score: (alpha0, alpha_X, alpha_V)

    pshat <- function(beta_hat,
                      data.obs) {
      # beta_hat <- betaEst_JICBPS
      # data.obs <-  unisimdataset

      beta_hat <- as.matrix(beta_hat)

      CovInfo_hat <- data.obs
      CovInfo_hat$Y <- NULL
      CovInfo_hat$W <- NULL
      CovInfo_hat <- as.matrix(cbind(1, CovInfo_hat))

      ps_hat <- 1 / (1 + exp(CovInfo_hat %*% (-beta_hat)))
      # There is probably an X'X missing somewhere that is causing these variance problems.
      probs.min <- 1e-6
      ps_hat <- pmin(1 - probs.min, ps_hat)
      ps_hat <- pmax(probs.min, ps_hat)
      ps_hat <- as.vector(ps_hat)

      return(ps_hat)
    }

    # ****** ------ correction propensity score ------ ****** #

    # The propensity score estimating function of DiShu
    # Input: the estimating coefficients of propensity score: (alpha0, alpha_X, alpha_V)

    pshat_DiShu <- function(beta_hat,
                            data.obs,
                            var_eps = var_eps) {

      # beta_hat <- betaEst_JICBPS
      # data.obs <-  unisimdataset
      # var_eps <- 0.1

      # para_beta
      para_X1_hat <- beta_hat["para_X1"]
      para_X2_hat <- beta_hat["para_X2"]
      para_X_hat <- rbind(para_X1_hat, para_X2_hat)

      CovInfo_hat <- data.obs
      CovInfo_hat$Y <- NULL
      CovInfo_hat$W <- NULL
      Xstar_hat <- CovInfo_hat[, c("X1", "X2")]
      V_hat <- CovInfo_hat[, c("V1", "V2")]

      # The measurement error covariance matrix containing only error-prone covariate.
      ifelse(dim_X == 1,
             Sigma_e <- var_eps,
             Sigma_e <- diag(x = rep(var_eps, dim_X)))

      delta_hat <- Xstar_hat +
        matrix(data = rep(Sigma_e %*% para_X_hat, each = N), nrow = N) * (W - 0.5)

      CovInfo_hat <- as.matrix(cbind(One = 1, delta = delta_hat, V = V_hat))
      ps_hat_star <- 1 / (1 + exp(CovInfo_hat %*% (-beta_hat)))
      # There is probably an X'X missing somewhere that is causing these variance problems.
      probs.min <- 1e-6
      ps_hat_star <- pmin(1 - probs.min, ps_hat_star)
      ps_hat_star <- pmax(probs.min, ps_hat_star)
      ps_hat_star <- as.vector(ps_hat_star)

      return(ps_hat_star)
    }

    # ****** ------ ATE with IPW ------ ****** #

    # The ATE estimating function for continuous Y
    # Input: the estimated propensity score

    ATE_IPW <- function(ps) {
      # ps <- ps_dfsane

      # Normalization
      Weights_ps_W1 <- W / ps
      Weights_ps_W0 <- (1 - W) / (1 - ps)
      Weights_ps <- Weights_ps_W1 / sum(Weights_ps_W1) + Weights_ps_W0 / sum(Weights_ps_W0)

      EY1_hat <- sum(Weights_ps * W * Y)
      EY0_hat <- sum(Weights_ps * (1 - W) * Y)
      estimate_ATE <- EY1_hat - EY0_hat

      # # paper version
      # EY1_star <- mean(W * Y / ps)
      # EY0_star <- mean((1 - W) * Y / (1 - ps))
      # estimate_ATE <- EY1_star - EY0_star

      return(estimate_ATE)
    }

    # ****** ------ mu with IPW ------ ****** #

    mu_IPW <- function(ps) {
      # The mu_IPW estimator follows CBPS paper in Page.251

      # ps <- ps_dfsane

      mu_hat <- sum(W * Y / ps) / sum(W / ps)

      return(mu_hat)
    }

    # ****** ------ Two-Step GMM ------ ****** #

    TwoStepGMM <- function(para_beta,
                           dataset,
                           objFunName = c("GMM.Sstar",
                                          "GMM.Huang",
                                          "GMM.Astar",
                                          "GMM.OIcCBPS")) {

      # objFunName <- "GMM.Sstar"
      # para_beta <- c(para_1 = 0,
      #                para_X1 = -1, para_X2 = 0.5,
      #                para_V1 = -0.25, para_V2 = -0.1)
      # dataset <- unisimdataset

      ## TwoStepGMM.loss - 写进函数里是为了方便并行计算
      TwoStepGMM.loss <- function(para_beta,
                                  dataset,
                                  invW.mat = NULL,
                                  objFunName) {

        # para_beta <- c(para_1 = 0,
        #                para_X1 = -1, para_X2 = 0.5,
        #                para_V1 = -0.25, para_V2 = -0.1)
        # dataset <- unisimdataset
        # objFunName <- "GMM.Sstar"

        g.mat <- eval(call(name = objFunName,
                           para_beta = para_beta,
                           x = dataset))
        gbar <- colMeans(g.mat)
        gbar <- as.matrix(gbar)

        if (is.null(invW.mat)) {
          invWhat_FirstStep <- diag(1, nrow = nrow(gbar), ncol = nrow(gbar))

          loss <- t(gbar) %*% invWhat_FirstStep %*% gbar
        }else {
          loss <- t(gbar) %*% invW.mat %*% gbar
        }

        return(loss)
      }

      opt.FirstStep <- optim.compiler(par = para_beta,
                                      fn = TwoStepGMM.loss,
                                      # control = list("maxit" = 1000, trace = T),
                                      control = list("maxit" = 1000),
                                      method = "BFGS",
                                      hessian = TRUE,
                                      dataset = dataset,
                                      objFunName = objFunName)
      # print(opt.FirstStep$par)

      thetahat.FirstStep <- opt.FirstStep$par

      g.mat.FirstStep <- eval(call(name = objFunName,
                                   para_beta = thetahat.FirstStep,
                                   x = dataset))
      g.mat.FirstStep <- t(g.mat.FirstStep)

      temp_What <- 0
      for (i in 1:N) {
        temp_What <- g.mat.FirstStep[, i] %*% t(g.mat.FirstStep[, i]) + temp_What
      }
      What_SecondStep <- temp_What / N

      # SVD_What_SecondStep <- svd(x = What_SecondStep)
      # SVD_What_SecondStep$d[SVD_What_SecondStep$d < 0] <- 0
      # inverse_What_SecondStep <- SVD_What_SecondStep$u %*%
      #   diag(1 / SVD_What_SecondStep$d) %*%
      #   t(SVD_What_SecondStep$v)
      # round(inverse_What_SecondStep, 3) == round(ginv(X = What_SecondStep), 3)
      invWhat_SecondStep <- ginv(X = What_SecondStep)

      opt.SecondStep <- optim.compiler(par = thetahat.FirstStep,
                                       fn = TwoStepGMM.loss,
                                       # control = list("maxit" = 1000, trace = T),
                                       control = list("maxit" = 1000),
                                       method = "BFGS",
                                       hessian = TRUE,
                                       dataset = dataset,
                                       invW.mat = invWhat_SecondStep,
                                       objFunName = objFunName)
      loss.SecondStep <- TwoStepGMM.loss(para_beta = thetahat.FirstStep,
                                         dataset = dataset,
                                         invW.mat = invWhat_SecondStep,
                                         objFunName = objFunName)

      thetahat.SecondStep <- opt.SecondStep$par
      ConvergeIndex <- (!opt.FirstStep$convergence) & (!opt.SecondStep$convergence)

      TwpStepGMM.result <- list(betaEst = thetahat.SecondStep,
                                ConvergeIndex = ConvergeIndex,
                                loss = loss.SecondStep,
                                invW = invWhat_SecondStep,
                                val = opt.SecondStep$value)

      return(TwpStepGMM.result)
    }

    # ****** ------ alpha.func ------ ****** #

    alpha.func <- function(alpha,
                           para_beta,
                           x,
                           objFunName) {

      # alpha <- 0.03
      # para_beta <- c(para_1 = 0,
      #                para_X1 = -1, para_X2 = 0.5,
      #                para_V1 = -0.25, para_V2 = -0.1)
      # x <- unisimdataset

      gmm.alpharesult <- TwoStepGMM(para_beta = alpha * para_beta,
                                    dataset = x,
                                    objFunName = objFunName)
      gmm.alphaloss <- gmm.alpharesult$loss

      return(gmm.alphaloss)
    }

    # ****** ------ ASMD 4 ks_data all predictors ------ ****** #

    # ASMD: ASMD Calculation Function

    ASMDCal4ksdata <- function(ps, # 1/Weights
                               Predictors_4Covs,
                               TrtIndex = TrtIndex, ColIndex = ColIndex) {

      # Predictors_4Covs <- Predictors_Covs_X
      # ps <- ps.mat[,"Real-ps"]

      # Normalization
      Weights_ps_W1 <- W / ps
      Weights_ps_W0 <- (1 - W) / (1 - ps)
      Weights <- Weights_ps_W1 / sum(Weights_ps_W1) + Weights_ps_W0 / sum(Weights_ps_W0)


      Z_NoME <- Predictors_4Covs
      Z1_NoME <- Predictors_4Covs[TrtIndex,]
      Z0_NoME <- Predictors_4Covs[ColIndex,]
      # Weighted Covariate (covariate without measurement error)
      # WeightedCovariate0 <- Z_NoME * matrix(rep(Weights, ncol(Z_NoME)), ncol = ncol(Z_NoME))
      WeightedCovariate <- Z_NoME * as.vector(Weights)

      # ASMD for each Covariate
      ## Numerator
      ASMD.Numerator_trt <- colSums(WeightedCovariate[TrtIndex,])
      ASMD.Numerator_col <- colSums(WeightedCovariate[ColIndex,])
      ## Denominator of ASMD 4 ATE
      ASMD.Denominator_partTrt <- apply(Z1_NoME, MARGIN = 2, sd)
      ASMD.Denominator_partCol <- apply(Z0_NoME, MARGIN = 2, sd)
      ASMD.Denominator <- sqrt((ASMD.Denominator_partTrt + ASMD.Denominator_partCol) / 2)

      ##ASMD
      ASMD <- abs(ASMD.Numerator_trt - ASMD.Numerator_col) / ASMD.Denominator

      return(data.frame(ASMD))
    }

  }

    # ====== Global Parameter Setting ====== #

  {

    # Name Variables
    MethodsName <- c("Real-ps",
                     "JICBPS", "OICBPS",
                     "DiShu",
                     "Huang",
                     "JIcCBPS", "OIcCBPS")
    paraName <- c("para_1", "para_X1", "para_X2", "para_V1", "para_V2")

    dim_X <- 2
    dim_V <- 2
    dim_Z <- dim_X + dim_V

    # There is probably an X'X missing somewhere that is causing these variance problems.
    probs.min <- 1e-6

  }

    Predictors <- data$Predictors
    Predictors_Covs_X <- Predictors[, c("x1", "x2", "x3", "x4")]
    Predictors_Covs_U <- Predictors[, c("u1", "u2", "u3", "u4")]
    unisimdataset <- as.data.frame(data$Dataset_obs)
    Y <- unisimdataset$Y
    W <- unisimdataset$W
    TrtIndex <- which(W == 1)
    ColIndex <- which(W == 0)
    Xstar <- unisimdataset[, c("X1", "X2")]
    V <- unisimdataset[, c("V1", "V2")]
    # The original propensity score of real data
    ps <- data$ps

    # var_eps <- varE
    # The measurement error covariance matrix containing only error-prone covariate.
    ifelse(dim_X == 1,
           Sigma_e <- var_eps,
           Sigma_e <- diag(x = rep(var_eps, dim_X)))
    # Sigma: the measurement error covariance matrix containing both error-prone and error-free covariate
    Sigma <- diag(c(rep(var_eps, dim_X), rep(0, dim_V)))

    # ========= naive estimate & Two-Step GMM Solutions ========= #

    formula4CBPS <- as.formula(W ~ X1 + X2 + V1 + V2)
    naive_Result <- CBPS(formula = formula4CBPS,
                         data = unisimdataset,
                         # Set to 0 to find the average treatment effect.
                         ATT = 0,
                         method = "exact",
                         verbose = FALSE)
    naive_beta <- as.vector(naive_Result$coefficients)
    names(naive_beta) <- paraName

    objFunName <- c("GMM.JICBPS", "GMM.OICBPS",
                    "GMM.Sstar", "GMM.Huang",
                    "GMM.Astar")

    TwoStepGMM.Result <- lapply(X = objFunName,
                                FUN = TwoStepGMM,
                                para_beta = naive_beta,
                                dataset = unisimdataset)
    names(TwoStepGMM.Result) <- objFunName

    # ========= naive JICBPS - just-identified ========= #

    betaEst_JICBPS <- TwoStepGMM.Result$GMM.JICBPS$betaEst

    # ========= naive OICBPS - over-identified ========= #

    betaEst_OICBPS <- TwoStepGMM.Result$GMM.OICBPS$betaEst

    # ========= DiShu ========= #

    betaEst_DiShu_GMM <- TwoStepGMM.Result$GMM.Sstar$betaEst

    # ========= Huang & Wang (2001) ========= #

    betaEst_Huang_GMM <- TwoStepGMM.Result$GMM.Huang$betaEst

    # ========= JIcCBPS ========= #

    betaEst_JIcCBPS <- TwoStepGMM.Result$GMM.Astar$betaEst

    # ========= OIcCBPS ========= #

    ##GLM estimation
    glm1 <- suppressWarnings(glm(formula = formula4CBPS,
                                 family = binomial,
                                 data = unisimdataset))
    glm1$coef[is.na(glm1$coef)] <- 0
    probs.glm <- glm1$fit
    glm1$fit <- probs.glm <- pmin(1 - probs.min, probs.glm)
    glm1$fit <- probs.glm <- pmax(probs.min, probs.glm)
    beta.curr <- glm1$coef
    beta.curr[is.na(beta.curr)] <- 0
    names(beta.curr) <- paraName

    # alpha.func <- function(alpha)gmm.loss(beta.curr * alpha)
    const <- optimize.compiler(alpha.func,
                               interval = c(.8, 1.1),
                               para_beta = beta.curr,
                               x = unisimdataset,
                               objFunName = "GMM.OIcCBPS")$minimum
    beta.curr <- beta.curr * const
    ##Generate estimates for balance and CBPSE
    gmm.init <- beta.curr

    gmm.glm.init <- TwoStepGMM(para_beta = gmm.init,
                               dataset = unisimdataset,
                               objFunName = "GMM.OIcCBPS")

    gmm.bal.init <- TwoStepGMM(para_beta = betaEst_JIcCBPS,
                               dataset = unisimdataset,
                               objFunName = "GMM.OIcCBPS")

    if (gmm.glm.init$val < gmm.bal.init$val) {
      GMMResults_OIcCBPS <- gmm.glm.init
    }else {
      GMMResults_OIcCBPS <- gmm.bal.init }

    betaEst_OIcCBPS <- GMMResults_OIcCBPS$betaEst

    # ------- betaEst ------ #

    beta.mat <- cbind(betaEst_JICBPS,
                      betaEst_OICBPS,
                      betaEst_DiShu_GMM,
                      betaEst_Huang_GMM,
                      betaEst_JIcCBPS,
                      betaEst_OIcCBPS)

    colnames(beta.mat) <- MethodsName[-1]
    rownames(beta.mat) <- paraName

    # ------- CvgIndex ------ #
    CvgIndex.mat <- cbind(TwoStepGMM.Result$GMM.JICBPS$ConvergeIndex,
                          TwoStepGMM.Result$GMM.OICBPS$ConvergeIndex,
                          TwoStepGMM.Result$GMM.Sstar$ConvergeIndex,
                          TwoStepGMM.Result$GMM.Huang$ConvergeIndex,
                          TwoStepGMM.Result$GMM.Astar$ConvergeIndex,
                          GMMResults_OIcCBPS$ConvergeIndex)
    colnames(CvgIndex.mat) <- MethodsName[-1]

    # ------- propensity score ------ #

    ps.mat.naiveCBPS <- apply(X = beta.mat[, c("JICBPS", "OICBPS")],
                              MARGIN = 2,
                              FUN = pshat,
                              data.obs = unisimdataset)

    ps.mat.otherMethods <- apply(X = beta.mat[, c("DiShu", "Huang", "JIcCBPS", "OIcCBPS")],
                                 MARGIN = 2,
                                 FUN = pshat_DiShu,
                                 data.obs = unisimdataset,
                                 var_eps = var_eps)

    ps.mat <- cbind(ps, ps.mat.naiveCBPS, ps.mat.otherMethods)
    colnames(ps.mat) <- MethodsName

    # ------- ASMD ------ #

    ## ASMD of predictors "x1, x2, x3, x4" from ks_data
    list_ASMD_Predictors_X <- apply(X = ps.mat, MARGIN = 2,
                                    FUN = ASMDCal4ksdata,
                                    Predictors_4Covs = Predictors_Covs_X,
                                    TrtIndex = TrtIndex, ColIndex = ColIndex)
    ASMD_Predictors_X <- matrix(data = unlist(x = list_ASMD_Predictors_X),
                                nrow = dim_Z,
                                byrow = F)
    colnames(ASMD_Predictors_X) <- MethodsName

    ## ASMD of predictors "u1, u2, u3, u4" from ks_data
    list_ASMD_Predictors_U <- apply(X = ps.mat, MARGIN = 2,
                                    FUN = ASMDCal4ksdata,
                                    Predictors_4Covs = Predictors_Covs_U,
                                    TrtIndex = TrtIndex, ColIndex = ColIndex)
    ASMD_Predictors_U <- matrix(data = unlist(x = list_ASMD_Predictors_U),
                                nrow = dim_Z,
                                byrow = F)
    colnames(ASMD_Predictors_U) <- MethodsName

    mat.ASMD <- rbind(ASMD_Predictors_X, ASMD_Predictors_U)
    rownames(mat.ASMD) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
    mat.ASMD <- as.data.frame(mat.ASMD)

    # ------- ATE ------ #

    ATE.mat <- apply(X = ps.mat,
                     MARGIN = 2,
                     FUN = ATE_IPW)

    # ------- Result ------ #

    unisim.ResultData <- list(unisim.beta = beta.mat,
                              unisim.CvgIndex = CvgIndex.mat,
                              unisim.ASMD = mat.ASMD,
                              unisim.ATE = ATE.mat)

    return(unisim.ResultData)
  }

}

}


# ====== Global Parameter Setting ====== #
{
  # Sample size
  N <- 1000
  # The number of simulations
  N_Sim <- 1000

  # Covariates Information
  dim_X <- 2
  dim_V <- 2
  dim_Z <- dim_X + dim_V

  # Variance of Measurement error
  # vec_sigmaE <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  # vec_varE <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  # vec_varE <- c(0.0, 0.1, 0.2, 0.3, 0.4)
  # vec_varE <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  # vec_varE <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  vec_varE <- c(0.0, 0.1, 0.2, 0.3, 0.5)

  # # use true covariates or not
  Y_scen <- "a"
  PS_scen <- "a"

  # Treatment and Propensity Score
  ## para_beta = (para_One, para_X, para_V)^T
  # True_beta <- as.matrix(c(para_1 = 0.2, para_X = 1, para_V = 1))
  # True_beta <- as.matrix(c(para_1 = 0.2, para_X = 0.5, para_V = 1))
  # True_beta <- as.matrix(c(para_1 = 0, para_X = -1, para_V = 0.5))
  True_beta <- as.matrix(c(para_1 = 0,
                           para_X1 = -1, para_X2 = 0.5,
                           para_V1 = -0.25, para_V2 = -0.1))

  # Outcome variable
  ## Covariance matrix of potential outcomes
  Sigma_Y <- 10
  ## Correlation coefficient matrix of potential outcomes
  rhoY01 <- 0
  ## True value of average treatment effect
  TrueValue_ATE <- 10
  ## True value of average effect of population
  TrueValue_mu <- 210

  AlgoIterMax <- 1000

  # There is probably an X'X missing somewhere that is causing these variance problems.
  probs.min <- 1e-6
}

# ========== Result Record Setting ========== #
{
  MethodsName <- c("Real-ps",
                   "JICBPS", "OICBPS",
                   "DiShu",
                   "Huang",
                   "JIcCBPS", "OIcCBPS")
  paraName <- c("para_1", "para_X1", "para_X2", "para_V1", "para_V2")
  varEName <- paste0("sigma_E^2: ", vec_varE)

  # ============ Result Table - beta ============ #

  # Mark down the beta estimating result
  EachSim_1000beta_JICBPS <- matrix(data = 0,
                                    nrow = N_Sim,
                                    ncol = dim_Z + 1)
  colnames(EachSim_1000beta_JICBPS) <- paraName
  EachSim_1000beta_OICBPS <- EachSim_1000beta_JICBPS
  EachSim_1000beta_DiShu <- EachSim_1000beta_JICBPS
  EachSim_1000beta_Huang <- EachSim_1000beta_JICBPS
  EachSim_1000beta_JIcCBPS <- EachSim_1000beta_JICBPS
  EachSim_1000beta_OIcCBPS <- EachSim_1000beta_JICBPS
  # EachSim_1000beta_Phi_MinusPlus_SumMinusPlus <- EachSim_1000beta_JICBPS

  # SCN: scenario
  ## just-identified CBPS
  AllSCNbetaEst_JICBPS <- list()
  AllSCNbetaEst_JICBPS$varE01 <- EachSim_1000beta_JICBPS
  AllSCNbetaEst_JICBPS$varE03 <- EachSim_1000beta_JICBPS
  AllSCNbetaEst_JICBPS$varE05 <- EachSim_1000beta_JICBPS
  AllSCNbetaEst_JICBPS$varE07 <- EachSim_1000beta_JICBPS
  AllSCNbetaEst_JICBPS$varE09 <- EachSim_1000beta_JICBPS
  ## over-identified CBPS
  AllSCNbetaEst_OICBPS <- list()
  AllSCNbetaEst_OICBPS$varE01 <- EachSim_1000beta_OICBPS
  AllSCNbetaEst_OICBPS$varE03 <- EachSim_1000beta_OICBPS
  AllSCNbetaEst_OICBPS$varE05 <- EachSim_1000beta_OICBPS
  AllSCNbetaEst_OICBPS$varE07 <- EachSim_1000beta_OICBPS
  AllSCNbetaEst_OICBPS$varE09 <- EachSim_1000beta_OICBPS
  ## DiShu
  AllSCNbetaEst_DiShu <- list()
  AllSCNbetaEst_DiShu$varE01 <- EachSim_1000beta_DiShu
  AllSCNbetaEst_DiShu$varE03 <- EachSim_1000beta_DiShu
  AllSCNbetaEst_DiShu$varE05 <- EachSim_1000beta_DiShu
  AllSCNbetaEst_DiShu$varE07 <- EachSim_1000beta_DiShu
  AllSCNbetaEst_DiShu$varE09 <- EachSim_1000beta_DiShu
  ## Huang
  AllSCNbetaEst_Huang <- list()
  AllSCNbetaEst_Huang$varE01 <- EachSim_1000beta_Huang
  AllSCNbetaEst_Huang$varE03 <- EachSim_1000beta_Huang
  AllSCNbetaEst_Huang$varE05 <- EachSim_1000beta_Huang
  AllSCNbetaEst_Huang$varE07 <- EachSim_1000beta_Huang
  AllSCNbetaEst_Huang$varE09 <- EachSim_1000beta_Huang
  ## JIcCBPS
  AllSCNbetaEst_JIcCBPS <- list()
  AllSCNbetaEst_JIcCBPS$varE01 <- EachSim_1000beta_JIcCBPS
  AllSCNbetaEst_JIcCBPS$varE03 <- EachSim_1000beta_JIcCBPS
  AllSCNbetaEst_JIcCBPS$varE05 <- EachSim_1000beta_JIcCBPS
  AllSCNbetaEst_JIcCBPS$varE07 <- EachSim_1000beta_JIcCBPS
  AllSCNbetaEst_JIcCBPS$varE09 <- EachSim_1000beta_JIcCBPS
  ## OIcCBPS
  AllSCNbetaEst_OIcCBPS <- list()
  AllSCNbetaEst_OIcCBPS$varE01 <- EachSim_1000beta_OIcCBPS
  AllSCNbetaEst_OIcCBPS$varE03 <- EachSim_1000beta_OIcCBPS
  AllSCNbetaEst_OIcCBPS$varE05 <- EachSim_1000beta_OIcCBPS
  AllSCNbetaEst_OIcCBPS$varE07 <- EachSim_1000beta_OIcCBPS
  AllSCNbetaEst_OIcCBPS$varE09 <- EachSim_1000beta_OIcCBPS

  # ============ Result Table - ASMD ============ #

  # cols of ASMD matrix: ASMD(predictors_X), ASMD(predictors_U)
  EachSim_ASMD_JICBPS <- matrix(data = 0,
                                ncol = dim_Z * 2, # ncol = x1 - x4,u1 - u4
                                nrow = N_Sim)
  colnames(EachSim_ASMD_JICBPS) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  EachSim_ASMD_RealPS <- EachSim_ASMD_JICBPS
  EachSim_ASMD_OICBPS <- EachSim_ASMD_JICBPS
  EachSim_ASMD_DiShu <- EachSim_ASMD_JICBPS
  EachSim_ASMD_Huang <- EachSim_ASMD_JICBPS
  EachSim_ASMD_JIcCBPS <- EachSim_ASMD_JICBPS
  EachSim_ASMD_OIcCBPS <- EachSim_ASMD_JICBPS
  # EachSim_ASMD_Phi_MinusPlus_SumMinusPlus <- EachSim_ASMD_JICBPS

  ## Real-PS
  AllSCNASMD_RealPS <- list()
  AllSCNASMD_RealPS$varE01 <- EachSim_ASMD_RealPS
  AllSCNASMD_RealPS$varE03 <- EachSim_ASMD_RealPS
  AllSCNASMD_RealPS$varE05 <- EachSim_ASMD_RealPS
  AllSCNASMD_RealPS$varE07 <- EachSim_ASMD_RealPS
  AllSCNASMD_RealPS$varE09 <- EachSim_ASMD_RealPS
  ## just-identified CBPS
  AllSCNASMD_JICBPS <- list()
  AllSCNASMD_JICBPS$varE01 <- EachSim_ASMD_JICBPS
  AllSCNASMD_JICBPS$varE03 <- EachSim_ASMD_JICBPS
  AllSCNASMD_JICBPS$varE05 <- EachSim_ASMD_JICBPS
  AllSCNASMD_JICBPS$varE07 <- EachSim_ASMD_JICBPS
  AllSCNASMD_JICBPS$varE09 <- EachSim_ASMD_JICBPS
  ## over-identified CBPS
  AllSCNASMD_OICBPS <- list()
  AllSCNASMD_OICBPS$varE01 <- EachSim_ASMD_OICBPS
  AllSCNASMD_OICBPS$varE03 <- EachSim_ASMD_OICBPS
  AllSCNASMD_OICBPS$varE05 <- EachSim_ASMD_OICBPS
  AllSCNASMD_OICBPS$varE07 <- EachSim_ASMD_OICBPS
  AllSCNASMD_OICBPS$varE09 <- EachSim_ASMD_OICBPS
  ## DiShu
  AllSCNASMD_DiShu <- list()
  AllSCNASMD_DiShu$varE01 <- EachSim_ASMD_DiShu
  AllSCNASMD_DiShu$varE03 <- EachSim_ASMD_DiShu
  AllSCNASMD_DiShu$varE05 <- EachSim_ASMD_DiShu
  AllSCNASMD_DiShu$varE07 <- EachSim_ASMD_DiShu
  AllSCNASMD_DiShu$varE09 <- EachSim_ASMD_DiShu
  ## Huang
  AllSCNASMD_Huang <- list()
  AllSCNASMD_Huang$varE01 <- EachSim_ASMD_Huang
  AllSCNASMD_Huang$varE03 <- EachSim_ASMD_Huang
  AllSCNASMD_Huang$varE05 <- EachSim_ASMD_Huang
  AllSCNASMD_Huang$varE07 <- EachSim_ASMD_Huang
  AllSCNASMD_Huang$varE09 <- EachSim_ASMD_Huang
  ## JIcCBPS
  AllSCNASMD_JIcCBPS <- list()
  AllSCNASMD_JIcCBPS$varE01 <- EachSim_ASMD_JIcCBPS
  AllSCNASMD_JIcCBPS$varE03 <- EachSim_ASMD_JIcCBPS
  AllSCNASMD_JIcCBPS$varE05 <- EachSim_ASMD_JIcCBPS
  AllSCNASMD_JIcCBPS$varE07 <- EachSim_ASMD_JIcCBPS
  AllSCNASMD_JIcCBPS$varE09 <- EachSim_ASMD_JIcCBPS
  ## OIcCBPS
  AllSCNASMD_OIcCBPS <- list()
  AllSCNASMD_OIcCBPS$varE01 <- EachSim_ASMD_OIcCBPS
  AllSCNASMD_OIcCBPS$varE03 <- EachSim_ASMD_OIcCBPS
  AllSCNASMD_OIcCBPS$varE05 <- EachSim_ASMD_OIcCBPS
  AllSCNASMD_OIcCBPS$varE07 <- EachSim_ASMD_OIcCBPS
  AllSCNASMD_OIcCBPS$varE09 <- EachSim_ASMD_OIcCBPS

  # ============ Result Table - ATE ============ #

  ## EachSim_1000ATE:
  # When we set the N_Sim = 1000, we run 1000 experiments under the same parameter setting.
  # We use EachSim_1000ATE to mark down the 1000 ATE estimates.
  EachSim_1000ATE <- matrix(data = 0,
                            nrow = N_Sim,
                            ncol = length(MethodsName))
  colnames(EachSim_1000ATE) <- MethodsName

  # Table for marking down all results in different scenarios
  AllSCNATEEst <- list()
  AllSCNATEEst$varE01 <- EachSim_1000ATE
  AllSCNATEEst$varE03 <- EachSim_1000ATE
  AllSCNATEEst$varE05 <- EachSim_1000ATE
  AllSCNATEEst$varE07 <- EachSim_1000ATE
  AllSCNATEEst$varE09 <- EachSim_1000ATE

  # ============ Result Table - Convergence ============ #

  # Mark the convergence
  UnCovergeMethodsName <- MethodsName[-1]
  EachSim_1000Cvg <- matrix(data = 0,
                            nrow = N_Sim,
                            ncol = length(UnCovergeMethodsName))
  colnames(EachSim_1000Cvg) <- UnCovergeMethodsName
  # Value = 1 indicates successful convergence.

  # Table for marking down all results in different scenarios
  AllSCNCvgEst <- list()
  AllSCNCvgEst$varE01 <- EachSim_1000Cvg
  AllSCNCvgEst$varE03 <- EachSim_1000Cvg
  AllSCNCvgEst$varE05 <- EachSim_1000Cvg
  AllSCNCvgEst$varE07 <- EachSim_1000Cvg
  AllSCNCvgEst$varE09 <- EachSim_1000Cvg

}

# ====== Simulation Setting ====== #

{
  No_TotalSim <- 0

  ith_varE <- 0

  # 计算可用线程数，并设置并行使用线程数
  no_cores <- detectCores() - 1

  # 初始化
  cl <- makeCluster(no_cores)
}

for (varE in vec_varE) {
  # varE <- 0.1
  # var_eps <- varE

  # Reset result table for each var_E scenario
{
  # Reset the matrix for marking down beta estimates
  EachSim_1000beta_JICBPS <- matrix(data = 0,
                                    nrow = N_Sim,
                                    ncol = dim_Z + 1)
  colnames(EachSim_1000beta_JICBPS) <- paraName
  EachSim_1000beta_OICBPS <- EachSim_1000beta_JICBPS
  EachSim_1000beta_DiShu <- EachSim_1000beta_JICBPS
  EachSim_1000beta_Huang <- EachSim_1000beta_JICBPS
  EachSim_1000beta_JIcCBPS <- EachSim_1000beta_JICBPS
  EachSim_1000beta_OIcCBPS <- EachSim_1000beta_JICBPS
  # EachSim_1000beta_Phi_MinusPlus_SumMinusPlus <- EachSim_1000beta_JICBPS

  # cols of ASMD matrix: ASMD(X), ASMD(V)
  # cols of ASMD matrix: ASMD(predictors_X), ASMD(predictors_U)
  EachSim_ASMD_JICBPS <- matrix(data = 0,
                                ncol = dim_Z * 2, # ncol = x1 - x4,u1 - u4
                                nrow = N_Sim)
  colnames(EachSim_ASMD_JICBPS) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  EachSim_ASMD_RealPS <- EachSim_ASMD_JICBPS
  EachSim_ASMD_OICBPS <- EachSim_ASMD_JICBPS
  EachSim_ASMD_DiShu <- EachSim_ASMD_JICBPS
  EachSim_ASMD_Huang <- EachSim_ASMD_JICBPS
  EachSim_ASMD_JIcCBPS <- EachSim_ASMD_JICBPS
  EachSim_ASMD_OIcCBPS <- EachSim_ASMD_JICBPS

  # Reset the matrix for marking down ATE estimates
  EachSim_1000ATE <- matrix(data = 0,
                            nrow = N_Sim,
                            ncol = length(MethodsName))
  colnames(EachSim_1000ATE) <- MethodsName


  # Reset matrix of the convergence rate
  EachSim_1000Cvg <- matrix(data = 0,
                            nrow = N_Sim,
                            ncol = length(UnCovergeMethodsName))
  colnames(EachSim_1000Cvg) <- UnCovergeMethodsName
}


  # Run N_sim simulation for specific var_E
{

  cat("\n # ================ # \n",
      "#  var_Eps = ", varE, " #",
      "\n # ================ # \n")

  AllSimData <- pblapply(cl = cl,
                         X = 1:N_Sim,
                         FUN = UniSimDataset,
                         N = N,
                         var_eps = varE)

  cat("\n # ============================== # \n",
      "  Dataset Construction Finished  ",
      "\n # ============================== # \n")

  # parallel::clusterExport(cl= cl,varlist = c("AllSimData"))

  AllSim.Result <- pblapply(cl = cl,
                            X = AllSimData,
                            FUN = UniSimATE,
                            N = N,
                            var_eps = varE)

  cat("\n # ========================== # \n",
      "  ATE Estimation Finished  ",
      "\n # ========================== # \n")

  for (ith_sim in 1:N_Sim) {

    # beta estimate
    EachSim_1000beta_JICBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.beta[, "JICBPS"]
    EachSim_1000beta_OICBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.beta[, "OICBPS"]
    EachSim_1000beta_DiShu[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.beta[, "DiShu"]
    EachSim_1000beta_Huang[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.beta[, "Huang"]
    EachSim_1000beta_JIcCBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.beta[, "JIcCBPS"]
    EachSim_1000beta_OIcCBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.beta[, "OIcCBPS"]

    # Convergence Index
    EachSim_1000Cvg[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.CvgIndex

    # ASMD
    ## Real-PS
    EachSim_ASMD_RealPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ASMD$`Real-ps`
    ## JICBPS
    EachSim_ASMD_JICBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ASMD$JICBPS
    ## OICBPS
    EachSim_ASMD_OICBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ASMD$OICBPS
    ## DiShu
    EachSim_ASMD_DiShu[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ASMD$DiShu
    ## Huang
    EachSim_ASMD_Huang[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ASMD$Huang
    ## JIcCBPS
    EachSim_ASMD_JIcCBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ASMD$JIcCBPS
    ## OIcCBPS
    EachSim_ASMD_OIcCBPS[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ASMD$OIcCBPS

    # ATE estimate
    EachSim_1000ATE[ith_sim,] <- AllSim.Result[[ith_sim]]$unisim.ATE

  }
}

  ith_varE <- ith_varE + 1

  if (ith_varE == 1) {
    AllSCNbetaEst_JICBPS[["varE01"]] <- EachSim_1000beta_JICBPS
    AllSCNbetaEst_OICBPS[["varE01"]] <- EachSim_1000beta_OICBPS
    AllSCNbetaEst_DiShu[["varE01"]] <- EachSim_1000beta_DiShu
    AllSCNbetaEst_Huang[["varE01"]] <- EachSim_1000beta_Huang
    AllSCNbetaEst_JIcCBPS[["varE01"]] <- EachSim_1000beta_JIcCBPS
    AllSCNbetaEst_OIcCBPS[["varE01"]] <- EachSim_1000beta_OIcCBPS

    AllSCNASMD_RealPS[["varE01"]] <- EachSim_ASMD_RealPS
    AllSCNASMD_JICBPS[["varE01"]] <- EachSim_ASMD_JICBPS
    AllSCNASMD_OICBPS[["varE01"]] <- EachSim_ASMD_OICBPS
    AllSCNASMD_DiShu[["varE01"]] <- EachSim_ASMD_DiShu
    AllSCNASMD_Huang[["varE01"]] <- EachSim_ASMD_Huang
    AllSCNASMD_JIcCBPS[["varE01"]] <- EachSim_ASMD_JIcCBPS
    AllSCNASMD_OIcCBPS[["varE01"]] <- EachSim_ASMD_OIcCBPS

    AllSCNATEEst[["varE01"]] <- EachSim_1000ATE
    AllSCNCvgEst[["varE01"]] <- EachSim_1000Cvg
  }else {
    if (ith_varE == 2) {
      AllSCNbetaEst_JICBPS[["varE03"]] <- EachSim_1000beta_JICBPS
      AllSCNbetaEst_OICBPS[["varE03"]] <- EachSim_1000beta_OICBPS
      AllSCNbetaEst_DiShu[["varE03"]] <- EachSim_1000beta_DiShu
      AllSCNbetaEst_Huang[["varE03"]] <- EachSim_1000beta_Huang
      AllSCNbetaEst_JIcCBPS[["varE03"]] <- EachSim_1000beta_JIcCBPS
      AllSCNbetaEst_OIcCBPS[["varE03"]] <- EachSim_1000beta_OIcCBPS

      AllSCNASMD_RealPS[["varE03"]] <- EachSim_ASMD_RealPS
      AllSCNASMD_JICBPS[["varE03"]] <- EachSim_ASMD_JICBPS
      AllSCNASMD_OICBPS[["varE03"]] <- EachSim_ASMD_OICBPS
      AllSCNASMD_DiShu[["varE03"]] <- EachSim_ASMD_DiShu
      AllSCNASMD_Huang[["varE03"]] <- EachSim_ASMD_Huang
      AllSCNASMD_JIcCBPS[["varE03"]] <- EachSim_ASMD_JIcCBPS
      AllSCNASMD_OIcCBPS[["varE03"]] <- EachSim_ASMD_OIcCBPS

      AllSCNATEEst[["varE03"]] <- EachSim_1000ATE
      AllSCNCvgEst[["varE03"]] <- EachSim_1000Cvg
    }else {
      if (ith_varE == 3) {
        AllSCNbetaEst_JICBPS[["varE05"]] <- EachSim_1000beta_JICBPS
        AllSCNbetaEst_OICBPS[["varE05"]] <- EachSim_1000beta_OICBPS
        AllSCNbetaEst_DiShu[["varE05"]] <- EachSim_1000beta_DiShu
        AllSCNbetaEst_Huang[["varE05"]] <- EachSim_1000beta_Huang
        AllSCNbetaEst_JIcCBPS[["varE05"]] <- EachSim_1000beta_JIcCBPS
        AllSCNbetaEst_OIcCBPS[["varE05"]] <- EachSim_1000beta_OIcCBPS

        AllSCNASMD_RealPS[["varE05"]] <- EachSim_ASMD_RealPS
        AllSCNASMD_JICBPS[["varE05"]] <- EachSim_ASMD_JICBPS
        AllSCNASMD_OICBPS[["varE05"]] <- EachSim_ASMD_OICBPS
        AllSCNASMD_DiShu[["varE05"]] <- EachSim_ASMD_DiShu
        AllSCNASMD_Huang[["varE05"]] <- EachSim_ASMD_Huang
        AllSCNASMD_JIcCBPS[["varE05"]] <- EachSim_ASMD_JIcCBPS
        AllSCNASMD_OIcCBPS[["varE05"]] <- EachSim_ASMD_OIcCBPS

        AllSCNATEEst[["varE05"]] <- EachSim_1000ATE
        AllSCNCvgEst[["varE05"]] <- EachSim_1000Cvg
      }else {
        if (ith_varE == 4) {
          AllSCNbetaEst_JICBPS[["varE07"]] <- EachSim_1000beta_JICBPS
          AllSCNbetaEst_OICBPS[["varE07"]] <- EachSim_1000beta_OICBPS
          AllSCNbetaEst_DiShu[["varE07"]] <- EachSim_1000beta_DiShu
          AllSCNbetaEst_Huang[["varE07"]] <- EachSim_1000beta_Huang
          AllSCNbetaEst_JIcCBPS[["varE07"]] <- EachSim_1000beta_JIcCBPS
          AllSCNbetaEst_OIcCBPS[["varE07"]] <- EachSim_1000beta_OIcCBPS

          AllSCNASMD_RealPS[["varE07"]] <- EachSim_ASMD_RealPS
          AllSCNASMD_JICBPS[["varE07"]] <- EachSim_ASMD_JICBPS
          AllSCNASMD_OICBPS[["varE07"]] <- EachSim_ASMD_OICBPS
          AllSCNASMD_DiShu[["varE07"]] <- EachSim_ASMD_DiShu
          AllSCNASMD_Huang[["varE07"]] <- EachSim_ASMD_Huang
          AllSCNASMD_JIcCBPS[["varE07"]] <- EachSim_ASMD_JIcCBPS
          AllSCNASMD_OIcCBPS[["varE07"]] <- EachSim_ASMD_OIcCBPS

          AllSCNATEEst[["varE07"]] <- EachSim_1000ATE
          AllSCNCvgEst[["varE07"]] <- EachSim_1000Cvg
        }else {
          AllSCNbetaEst_JICBPS[["varE09"]] <- EachSim_1000beta_JICBPS
          AllSCNbetaEst_OICBPS[["varE09"]] <- EachSim_1000beta_OICBPS
          AllSCNbetaEst_DiShu[["varE09"]] <- EachSim_1000beta_DiShu
          AllSCNbetaEst_Huang[["varE09"]] <- EachSim_1000beta_Huang
          AllSCNbetaEst_JIcCBPS[["varE09"]] <- EachSim_1000beta_JIcCBPS
          AllSCNbetaEst_OIcCBPS[["varE09"]] <- EachSim_1000beta_OIcCBPS

          AllSCNASMD_RealPS[["varE09"]] <- EachSim_ASMD_RealPS
          AllSCNASMD_JICBPS[["varE09"]] <- EachSim_ASMD_JICBPS
          AllSCNASMD_OICBPS[["varE09"]] <- EachSim_ASMD_OICBPS
          AllSCNASMD_DiShu[["varE09"]] <- EachSim_ASMD_DiShu
          AllSCNASMD_Huang[["varE09"]] <- EachSim_ASMD_Huang
          AllSCNASMD_JIcCBPS[["varE09"]] <- EachSim_ASMD_JIcCBPS
          AllSCNASMD_OIcCBPS[["varE09"]] <- EachSim_ASMD_OIcCBPS

          AllSCNATEEst[["varE09"]] <- EachSim_1000ATE
          AllSCNCvgEst[["varE09"]] <- EachSim_1000Cvg
        }
      }
    }
  }


}


# ------------- Result ------------- #

# ========== beta estimate ========== #

## JICBPS
{
  avgbeta_JICBPS <- matrix(unlist(lapply(X = AllSCNbetaEst_JICBPS,
                                         FUN = colMeans, na.rm = T)),
                           nrow = dim_Z + 1)
  bias_beta_JICBPS <- sweep(x = avgbeta_JICBPS,
                            MARGIN = 1,
                            STATS = True_beta,
                            FUN = "-")

  ReBias_beta_JICBPS <- apply(X = bias_beta_JICBPS, MARGIN = 2,
                              FUN = function(x) { x / True_beta * 100 })

  RMSE_beta_JICBPS <- matrix(data = unlist(lapply(X = AllSCNbetaEst_JICBPS,
                                                  FUN = function(x) {
                                                    # x <- AllSCNbetaEst_JICBPS[[1]]
                                                    RMSE <- NULL
                                                    for (i in 1:length(paraName)) {
                                                      RMSE[i] <- RMSEfun(actual = True_beta[i], predicted = x[, i])
                                                    }
                                                    return(RMSE) })),
                             nrow = dim_Z + 1)

  avgbeta_JICBPS <- round(x = avgbeta_JICBPS, digits = 5)

  AllSCN_avgbeta_JICBPS <- cbind(True_beta, avgbeta_JICBPS)
  rownames(AllSCN_avgbeta_JICBPS) <- paraName
  rownames(bias_beta_JICBPS) <- paraName
  rownames(ReBias_beta_JICBPS) <- paraName
  rownames(RMSE_beta_JICBPS) <- paraName

  colnames(AllSCN_avgbeta_JICBPS) <- c("TrueValue",
                                       # paste0("sigmaE", c("01", "03", "05", "07", "09"))
                                       varEName)
  colnames(bias_beta_JICBPS) <- varEName
  colnames(ReBias_beta_JICBPS) <- varEName
  colnames(RMSE_beta_JICBPS) <- varEName
}

## OICBPS
{
  avgbeta_OICBPS <- matrix(unlist(lapply(X = AllSCNbetaEst_OICBPS,
                                         FUN = colMeans, na.rm = T)),
                           nrow = dim_Z + 1)
  bias_beta_OICBPS <- sweep(x = avgbeta_OICBPS,
                            MARGIN = 1,
                            STATS = True_beta,
                            FUN = "-")

  ReBias_beta_OICBPS <- apply(X = bias_beta_OICBPS, MARGIN = 2,
                              FUN = function(x) { x / True_beta * 100 })

  RMSE_beta_OICBPS <- matrix(data = unlist(lapply(X = AllSCNbetaEst_OICBPS,
                                                  FUN = function(x) {
                                                    RMSE <- NULL
                                                    for (i in 1:length(paraName)) {
                                                      RMSE[i] <- RMSEfun(actual = True_beta[i], predicted = x[, i])
                                                    }
                                                    return(RMSE) })),
                             nrow = dim_Z + 1)

  avgbeta_OICBPS <- round(x = avgbeta_OICBPS, digits = 5)

  AllSCN_avgbeta_OICBPS <- cbind(True_beta, avgbeta_OICBPS)
  rownames(AllSCN_avgbeta_OICBPS) <- paraName
  rownames(bias_beta_OICBPS) <- paraName
  rownames(ReBias_beta_OICBPS) <- paraName
  rownames(RMSE_beta_OICBPS) <- paraName

  colnames(AllSCN_avgbeta_OICBPS) <- c("TrueValue",
                                       # paste0("sigmaE", c("01", "03", "05", "07", "09"))
                                       varEName)
  colnames(bias_beta_OICBPS) <- varEName
  colnames(ReBias_beta_OICBPS) <- varEName
  colnames(RMSE_beta_OICBPS) <- varEName
}

## DiShu
{
  # avgbeta_DiShu <- lapply(X = AllSCNbetaEst_DiShu,
  #                         FUN = colMeans, na.rm = T)
  avgbeta_DiShu <- matrix(unlist(lapply(X = AllSCNbetaEst_DiShu,
                                        FUN = colMeans, na.rm = T)),
                          nrow = dim_Z + 1)

  bias_beta_DiShu <- sweep(x = avgbeta_DiShu,
                           MARGIN = 1,
                           STATS = True_beta,
                           FUN = "-")

  ReBias_beta_DiShu <- apply(X = bias_beta_DiShu,
                             MARGIN = 2,
                             FUN = function(x) { x / True_beta * 100 })

  RMSE_beta_DiShu <- matrix(data = unlist(lapply(X = AllSCNbetaEst_DiShu,
                                                 FUN = function(x) {
                                                   RMSE <- NULL
                                                   for (i in 1:length(paraName)) {
                                                     RMSE[i] <- RMSEfun(actual = True_beta[i], predicted = x[, i])
                                                   }
                                                   return(RMSE) })),
                            nrow = dim_Z + 1)

  avgbeta_DiShu <- round(x = avgbeta_DiShu, digits = 5)

  AllSCN_avgbeta_DiShu <- cbind(True_beta, avgbeta_DiShu)
  rownames(AllSCN_avgbeta_DiShu) <- paraName
  rownames(bias_beta_DiShu) <- paraName
  rownames(ReBias_beta_DiShu) <- paraName
  rownames(RMSE_beta_DiShu) <- paraName

  colnames(AllSCN_avgbeta_DiShu) <- c("TrueValue",
                                      # paste0("sigmaE", c("01", "03", "05", "07", "09"))
                                      varEName)
  colnames(bias_beta_DiShu) <- varEName
  colnames(ReBias_beta_DiShu) <- varEName
  colnames(RMSE_beta_DiShu) <- varEName
}

## Huang
{
  # avgbeta_Huang <- lapply(X = AllSCNbetaEst_Huang,
  #                        FUN = colMeans, na.rm = T)
  avgbeta_Huang <- matrix(unlist(lapply(X = AllSCNbetaEst_Huang,
                                        FUN = colMeans, na.rm = T)),
                          nrow = dim_Z + 1)
  bias_beta_Huang <- sweep(x = avgbeta_Huang,
                           MARGIN = 1,
                           STATS = True_beta,
                           FUN = "-")

  ReBias_beta_Huang <- apply(X = bias_beta_Huang, MARGIN = 2,
                             FUN = function(x) { x / True_beta * 100 })

  RMSE_beta_Huang <- matrix(data = unlist(lapply(X = AllSCNbetaEst_Huang,
                                                 FUN = function(x) {
                                                   RMSE <- NULL
                                                   for (i in 1:length(paraName)) {
                                                     RMSE[i] <- RMSEfun(actual = True_beta[i], predicted = x[, i])
                                                   }
                                                   return(RMSE) })),
                            nrow = dim_Z + 1)

  avgbeta_Huang <- round(x = avgbeta_Huang, digits = 5)

  AllSCN_avgbeta_Huang <- cbind(True_beta, avgbeta_Huang)
  rownames(AllSCN_avgbeta_Huang) <- paraName
  rownames(bias_beta_Huang) <- paraName
  rownames(ReBias_beta_Huang) <- paraName
  rownames(RMSE_beta_Huang) <- paraName

  colnames(AllSCN_avgbeta_Huang) <- c("TrueValue",
                                      # paste0("sigmaE", c("01", "03", "05", "07", "09"))
                                      varEName)
  colnames(bias_beta_Huang) <- varEName
  colnames(ReBias_beta_Huang) <- varEName
  colnames(RMSE_beta_Huang) <- varEName
}

## JIcCBPS
{
  # avgbeta_JIcCBPS <- lapply(X = AllSCNbetaEst_JIcCBPS,
  #                        FUN = colMeans, na.rm = T)
  avgbeta_JIcCBPS <- matrix(unlist(lapply(X = AllSCNbetaEst_JIcCBPS,
                                          FUN = colMeans, na.rm = T)),
                            nrow = dim_Z + 1)
  bias_beta_JIcCBPS <- sweep(x = avgbeta_JIcCBPS,
                             MARGIN = 1,
                             STATS = True_beta,
                             FUN = "-")

  ReBias_beta_JIcCBPS <- apply(X = bias_beta_JIcCBPS, MARGIN = 2,
                               FUN = function(x) { x / True_beta * 100 })

  RMSE_beta_JIcCBPS <- matrix(data = unlist(lapply(X = AllSCNbetaEst_JIcCBPS,
                                                   FUN = function(x) {
                                                     RMSE <- NULL
                                                     for (i in 1:length(paraName)) {
                                                       RMSE[i] <- RMSEfun(actual = True_beta[i], predicted = x[, i])
                                                     }
                                                     return(RMSE) })),
                              nrow = dim_Z + 1)

  avgbeta_JIcCBPS <- round(x = avgbeta_JIcCBPS, digits = 5)

  AllSCN_avgbeta_JIcCBPS <- cbind(True_beta, avgbeta_JIcCBPS)
  rownames(AllSCN_avgbeta_JIcCBPS) <- paraName
  rownames(bias_beta_JIcCBPS) <- paraName
  rownames(ReBias_beta_JIcCBPS) <- paraName
  rownames(RMSE_beta_JIcCBPS) <- paraName

  colnames(AllSCN_avgbeta_JIcCBPS) <- c("TrueValue",
                                        # paste0("sigmaE", c("01", "03", "05", "07", "09"))
                                        varEName)
  colnames(bias_beta_JIcCBPS) <- varEName
  colnames(ReBias_beta_JIcCBPS) <- varEName
  colnames(RMSE_beta_JIcCBPS) <- varEName
}

## OIcCBPS
{
  # avgbeta_OIcCBPS <- lapply(X = AllSCNbetaEst_OIcCBPS,
  #                        FUN = colMeans, na.rm = T)
  avgbeta_OIcCBPS <- matrix(unlist(lapply(X = AllSCNbetaEst_OIcCBPS,
                                          FUN = colMeans, na.rm = T)),
                            nrow = dim_Z + 1)
  bias_beta_OIcCBPS <- sweep(x = avgbeta_OIcCBPS,
                             MARGIN = 1,
                             STATS = True_beta,
                             FUN = "-")

  ReBias_beta_OIcCBPS <- apply(X = bias_beta_OIcCBPS, MARGIN = 2,
                               FUN = function(x) { x / True_beta * 100 })

  RMSE_beta_OIcCBPS <- matrix(data = unlist(lapply(X = AllSCNbetaEst_OIcCBPS,
                                                   FUN = function(x) {
                                                     RMSE <- NULL
                                                     for (i in 1:length(paraName)) {
                                                       RMSE[i] <- RMSEfun(actual = True_beta[i], predicted = x[, i])
                                                     }
                                                     return(RMSE) })),
                              nrow = dim_Z + 1)

  avgbeta_OIcCBPS <- round(x = avgbeta_OIcCBPS, digits = 5)

  AllSCN_avgbeta_OIcCBPS <- cbind(True_beta, avgbeta_OIcCBPS)
  rownames(AllSCN_avgbeta_OIcCBPS) <- paraName
  rownames(bias_beta_OIcCBPS) <- paraName
  rownames(ReBias_beta_OIcCBPS) <- paraName
  rownames(RMSE_beta_OIcCBPS) <- paraName

  colnames(AllSCN_avgbeta_OIcCBPS) <- c("TrueValue",
                                        # paste0("sigmaE", c("01", "03", "05", "07", "09"))
                                        varEName)
  colnames(bias_beta_OIcCBPS) <- varEName
  colnames(ReBias_beta_OIcCBPS) <- varEName
  colnames(RMSE_beta_OIcCBPS) <- varEName
}

# ========== ASMD estimate ========== #

{
  ## RealPS
  avgASMD_RealPS <- t(matrix(unlist(lapply(X = AllSCNASMD_RealPS,
                                           FUN = colMeans, na.rm = T)),
                             nrow = dim_Z * 2))
  avgASMD_RealPS <- round(x = avgASMD_RealPS, digits = 3)
  rownames(avgASMD_RealPS) <- varEName
  colnames(avgASMD_RealPS) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  avgASMD_RealPS <- data.frame(sigma_epsilon2 = vec_varE, Method = "Real-ps", avgASMD_RealPS)

  ## just-identified CBPS
  avgASMD_JICBPS <- t(matrix(unlist(lapply(X = AllSCNASMD_JICBPS,
                                           FUN = colMeans, na.rm = T)),
                             nrow = dim_Z * 2))
  avgASMD_JICBPS <- round(x = avgASMD_JICBPS, digits = 3)
  rownames(avgASMD_JICBPS) <- varEName
  colnames(avgASMD_JICBPS) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  avgASMD_JICBPS <- data.frame(sigma_epsilon2 = vec_varE, Method = "JICBPS", avgASMD_JICBPS)

  ## over-identified CBPS
  avgASMD_OICBPS <- t(matrix(unlist(lapply(X = AllSCNASMD_OICBPS,
                                           FUN = colMeans, na.rm = T)),
                             nrow = dim_Z * 2))
  avgASMD_OICBPS <- round(x = avgASMD_OICBPS, digits = 3)
  rownames(avgASMD_OICBPS) <- varEName
  colnames(avgASMD_OICBPS) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  avgASMD_OICBPS <- data.frame(sigma_epsilon2 = vec_varE, Method = "OICBPS", avgASMD_OICBPS)

  ## DiShu
  avgASMD_DiShu <- t(matrix(unlist(lapply(X = AllSCNASMD_DiShu,
                                          FUN = colMeans, na.rm = T)),
                            nrow = dim_Z * 2))
  avgASMD_DiShu <- round(x = avgASMD_DiShu, digits = 3)
  rownames(avgASMD_DiShu) <- varEName
  colnames(avgASMD_DiShu) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  avgASMD_DiShu <- data.frame(sigma_epsilon2 = vec_varE, Method = "DiShu", avgASMD_DiShu)

  ## Huang
  avgASMD_Huang <- t(matrix(unlist(lapply(X = AllSCNASMD_Huang,
                                          FUN = colMeans, na.rm = T)),
                            nrow = dim_Z * 2))
  avgASMD_Huang <- round(x = avgASMD_Huang, digits = 3)
  rownames(avgASMD_Huang) <- varEName
  colnames(avgASMD_Huang) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  avgASMD_Huang <- data.frame(sigma_epsilon2 = vec_varE, Method = "Huang", avgASMD_Huang)

  ## JIcCBPS
  avgASMD_JIcCBPS <- t(matrix(unlist(lapply(X = AllSCNASMD_JIcCBPS,
                                            FUN = colMeans, na.rm = T)),
                              nrow = dim_Z * 2))
  avgASMD_JIcCBPS <- round(x = avgASMD_JIcCBPS, digits = 3)
  rownames(avgASMD_JIcCBPS) <- varEName
  colnames(avgASMD_JIcCBPS) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  avgASMD_JIcCBPS <- data.frame(sigma_epsilon2 = vec_varE, Method = "JIcCBPS", avgASMD_JIcCBPS)

  ## OIcCBPS
  avgASMD_OIcCBPS <- t(matrix(unlist(lapply(X = AllSCNASMD_OIcCBPS,
                                            FUN = colMeans, na.rm = T)),
                              nrow = dim_Z * 2))
  avgASMD_OIcCBPS <- round(x = avgASMD_OIcCBPS, digits = 3)
  rownames(avgASMD_OIcCBPS) <- varEName
  colnames(avgASMD_OIcCBPS) <- paste0(rep(paste0("ASMD_", c("x", "u")), each = dim_Z), 1:dim_Z)
  avgASMD_OIcCBPS <- data.frame(sigma_epsilon2 = vec_varE, Method = "OIcCBPS", avgASMD_OIcCBPS)


  mat_avgASMD <- rbind(avgASMD_RealPS,
                       avgASMD_JICBPS, avgASMD_OICBPS,
                       avgASMD_DiShu, avgASMD_Huang,
                       avgASMD_JIcCBPS, avgASMD_OIcCBPS)
}

# ========== ATE estimate ========== #

{
  avgATE <- matrix(data = unlist(lapply(X = AllSCNATEEst,
                                        FUN = colMeans,
                                        na.rm = T)),
                   ncol = length(MethodsName),
                   byrow = T)
  colnames(avgATE) <- MethodsName
  rownames(avgATE) <- varEName

  # Bias
  bias_ATE <- avgATE - TrueValue_ATE
  colnames(bias_ATE) <- MethodsName

  # ReBias%: the average relative bias in percent
  ReBias_ATE <- bias_ATE / TrueValue_ATE * 100
  colnames(ReBias_ATE) <- MethodsName

  # RMSE: root mean square error
  RMSE_ATE <- matrix(data = unlist(lapply(X = AllSCNATEEst,
                                          FUN = function(x) {
                                            apply(X = x, MARGIN = 2,
                                                  FUN = RMSEfun,
                                                  actual = TrueValue_ATE) })),
                     ncol = length(MethodsName),
                     byrow = T)
  colnames(RMSE_ATE) <- MethodsName
  rownames(RMSE_ATE) <- varEName

  # Real-ps
  Result_ATE_realps <- cbind(avgATE[, "Real-ps"],
                             bias_ATE[, "Real-ps"],
                             ReBias_ATE[, "Real-ps"],
                             RMSE_ATE[, "Real-ps"])
  colnames(Result_ATE_realps) <- c("ATE", "Bias", "ReBias", "RMSE")

  # JICBPS
  Result_ATE_JICBPS <- cbind(avgATE[, "JICBPS"],
                             bias_ATE[, "JICBPS"],
                             ReBias_ATE[, "JICBPS"],
                             RMSE_ATE[, "JICBPS"])
  colnames(Result_ATE_JICBPS) <- c("ATE", "Bias", "ReBias", "RMSE")

  # OICBPS
  Result_ATE_OICBPS <- cbind(avgATE[, "OICBPS"],
                             bias_ATE[, "OICBPS"],
                             ReBias_ATE[, "OICBPS"],
                             RMSE_ATE[, "OICBPS"])
  colnames(Result_ATE_OICBPS) <- c("ATE", "Bias", "ReBias", "RMSE")

  # DiShu
  Result_ATE_DiShu <- cbind(avgATE[, "DiShu"],
                            bias_ATE[, "DiShu"],
                            ReBias_ATE[, "DiShu"],
                            RMSE_ATE[, "DiShu"])
  colnames(Result_ATE_DiShu) <- c("ATE", "Bias", "ReBias", "RMSE")


  # Huang
  Result_ATE_Huang <- cbind(avgATE[, "Huang"],
                            bias_ATE[, "Huang"],
                            ReBias_ATE[, "Huang"],
                            RMSE_ATE[, "Huang"])
  colnames(Result_ATE_Huang) <- c("ATE", "Bias", "ReBias", "RMSE")

  # JIcCBPS
  Result_ATE_JIcCBPS <- cbind(avgATE[, "JIcCBPS"],
                              bias_ATE[, "JIcCBPS"],
                              ReBias_ATE[, "JIcCBPS"],
                              RMSE_ATE[, "JIcCBPS"])
  colnames(Result_ATE_JIcCBPS) <- c("ATE", "Bias", "ReBias", "RMSE")

  # OIcCBPS
  Result_ATE_OIcCBPS <- cbind(avgATE[, "OIcCBPS"],
                              bias_ATE[, "OIcCBPS"],
                              ReBias_ATE[, "OIcCBPS"],
                              RMSE_ATE[, "OIcCBPS"])
  colnames(Result_ATE_OIcCBPS) <- c("ATE", "Bias", "ReBias", "RMSE")
}

# ========== Convergence estimate ========== #

{

  AllSCNCvgRate <- matrix(unlist(lapply(X = AllSCNCvgEst, FUN = colMeans)),
                          ncol = length(UnCovergeMethodsName),
                          byrow = T)
  colnames(AllSCNCvgRate) <- UnCovergeMethodsName
  rownames(AllSCNCvgRate) <- varEName
}

# ========== Print beta Result ========== #
# ========== Print beta Result ========== #
# ========== Print beta Result ========== #
# ========== Print beta Result ========== #
# ========== Print beta Result ========== #

{
  ## just-identified CBPS
  cat("\n # ========== Print beta Result - JICBPS ========== #
   The average estimated beta: \n")
  print(round(AllSCN_avgbeta_JICBPS, digits = 3))

  cat("\n The bias of estimated beta: \n")
  print(round(bias_beta_JICBPS, digits = 3))

  cat("\n The average relative bias (ReBias) of beta in percent (%): \n")
  print(round(ReBias_beta_JICBPS, digits = 3))

  cat("\n The root mean square error (RMSE) of beta: \n")
  print(round(RMSE_beta_JICBPS, digits = 3))

  ## over-identified CBPS
  cat("\n # ========== Print beta Result - OICBPS ========== #
   The average estimated beta: \n")
  print(round(AllSCN_avgbeta_OICBPS, digits = 3))

  cat("\n The bias of estimated beta: \n")
  print(round(bias_beta_OICBPS, digits = 3))

  cat("\n The average relative bias (ReBias) of beta in percent (%): \n")
  print(round(ReBias_beta_OICBPS, digits = 3))

  cat("\n The root mean square error (RMSE) of beta: \n")
  print(round(RMSE_beta_OICBPS, digits = 3))

  ## DiShu
  cat("\n # ========== Print beta Result - DiShu ========== #
   The average estimated beta: \n")
  print(round(AllSCN_avgbeta_DiShu, digits = 3))

  cat("\n The bias of estimated beta: \n")
  print(round(bias_beta_DiShu, digits = 3))

  cat("\n The average relative bias (ReBias) of beta in percent (%): \n")
  print(round(ReBias_beta_DiShu, digits = 3))

  cat("\n The root mean square error (RMSE) of beta: \n")
  print(round(RMSE_beta_DiShu, digits = 3))

  ## Huang
  cat("\n # ========== Print beta Result - Huang ========== #
   The average estimated beta: \n")
  print(round(AllSCN_avgbeta_Huang, digits = 3))

  cat("\n The bias of estimated beta: \n")
  print(round(bias_beta_Huang, digits = 3))

  cat("\n The average relative bias (ReBias) of beta in percent (%): \n")
  print(round(ReBias_beta_Huang, digits = 3))

  cat("\n The root mean square error (RMSE) of beta: \n")
  print(round(RMSE_beta_Huang, digits = 3))

  ## JIcCBPS
  cat("\n # ========== Print beta Result - JIcCBPS ========== #
   The average estimated beta: \n")
  print(round(AllSCN_avgbeta_JIcCBPS, digits = 3))

  cat("\n The bias of estimated beta: \n")
  print(round(bias_beta_JIcCBPS, digits = 3))

  cat("\n The average relative bias (ReBias) of beta in percent (%): \n")
  print(round(ReBias_beta_JIcCBPS, digits = 3))

  cat("\n The root mean square error (RMSE) of beta: \n")
  print(round(RMSE_beta_JIcCBPS, digits = 3))

  ## OIcCBPS
  cat("\n # ========== Print beta Result - OIcCBPS ========== #
   The average estimated beta: \n")
  print(round(AllSCN_avgbeta_OIcCBPS, digits = 3))

  cat("\n The bias of estimated beta: \n")
  print(round(bias_beta_OIcCBPS, digits = 3))

  cat("\n The average relative bias (ReBias) of beta in percent (%): \n")
  print(round(ReBias_beta_OIcCBPS, digits = 3))

  cat("\n The root mean square error (RMSE) of beta: \n")
  print(round(RMSE_beta_OIcCBPS, digits = 3))

}

# ========== Print ATE Result ========== #
{
  cat("\n # === The ATE estimated with the true propensity score  === # \n")
  print(round(Result_ATE_realps, digits = 3))

  cat("\n # === The ATE estimated with the naive JICBPS  === # \n")
  print(round(Result_ATE_JICBPS, digits = 3))

  cat("\n # === The ATE estimated with the naive OICBPS  === # \n")
  print(round(Result_ATE_OICBPS, digits = 3))

  cat("\n # === The ATE estimated with DiShu  === # \n")
  print(round(Result_ATE_DiShu, digits = 3))

  cat("\n # === The ATE estimated with Huang  === # \n")
  print(round(Result_ATE_Huang, digits = 3))

  cat("\n # === The ATE estimated with JIcCBPS  === # \n")
  print(round(Result_ATE_JIcCBPS, digits = 3))

  cat("\n # === The ATE estimated with OIcCBPS  === # \n")
  print(round(Result_ATE_OIcCBPS, digits = 3))

  ATE_ALlMethod <- round(x = rbind(Result_ATE_realps,
                                   Result_ATE_JICBPS,
                                   Result_ATE_OICBPS,
                                   Result_ATE_DiShu,
                                   Result_ATE_Huang,
                                   Result_ATE_JIcCBPS,
                                   Result_ATE_OIcCBPS),
                         digits = 3)
  # Method <- rep(x = MethodsName, each = length(vec_sigmaE))
  Method <- rep(x = MethodsName, each = length(vec_varE))
  sigma_epsilon2 <- rep(x = vec_varE, length(MethodsName))
  ATE_ALlMethods <- data.frame(sigma_epsilon2, Method, ATE_ALlMethod,
                               row.names = NULL)
  ATE_ALlMethods <- ATE_ALlMethods[order(ATE_ALlMethods$sigma_epsilon2),]

  ATE_ALlMethods <- merge(x = ATE_ALlMethods, y = mat_avgASMD,
                          by.x = c("sigma_epsilon2", "Method"))

  csvfileName <- paste0("KSdata-4Covs_ATE_Y-a_PS-a_N", N, "_Sim", N_Sim, "_WeightsSum1.csv")

  write.csv(x = ATE_ALlMethods,
            file = csvfileName,
            row.names = FALSE)

  cat("\n # === The Convergence Rate of ATE estimates  === # \n")
  print(round(AllSCNCvgRate, digits = 3))
}

stopCluster(cl)

Time.endpoint <- Sys.time()

Time.elapse <- Time.endpoint - Time.startpoint

print(Time.elapse)

save.image("17.2.KSdata-4Covs_ATE_Y-a_PS-a.Rdata")