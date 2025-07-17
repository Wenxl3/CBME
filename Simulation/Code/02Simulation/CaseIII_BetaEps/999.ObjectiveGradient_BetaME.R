# -*- coding: utf-8 -*-
# ------------------------------------------
# Project     : EBME_R_202306022228
# Title       : 999.ObjectiveGradient_BetaME.R
# Objective   : Objective, gradient, hessian functions for mixture normal ME
# Reference   : (1) Yijian Huang & C. Y. Wang: Cox Regression with Accurate Covariates Unascertainable: A Nonparametric-Correction Approach
#               (2) Chengcheng Hu & D. Y Lin: Semiparametric Failure Time Regression With Replicates of Mismeasured Covariates
# Remark      : (1) Only the estimating equation for mismeasured covariates are constructed with the method of Huang & Wang.
#               (2) The estimating equation for precisely measuered covariates are constructed with the method of Hu & Lin.
#               (3) All units have the same number of replicats, which is set to be r = 2.
#               (4) CEB and BCEB use and only use the firsr replicate dataset.
#               (4) In simulation, CEB and BCEB use [TRUE] parameters in therir estimating equations, including the variance-covariance matrix $\Sigma$ of measurement error.
# Created on  : 2023/06/02 - 22:28
# ------------------------------------------

# ======================================== #
#           Method of naive EB             #
# ======================================== #

# EB: Entropy Balancing (Jens Hainmueller)----
objective.zqy.r_Beta <- function(theta) {
  f <- log(sum(drop((1 - W) * exp(Ze_Beta[[1]] %*% theta)))) +
    sum(-Ze1_bar_Beta[[1]] * theta)
  return(f)
}

gradient.zqy.r_Beta <- function(theta) {
  w <- drop((1 - W) * exp(Ze_Beta[[1]] %*% theta))
  f <- drop(w %*% (Ze_Beta[[1]]) / sum(w) - Ze1_bar_Beta[[1]])
  return(f)
}

hessian.zqy.r_Beta <- function(theta) {
  f <- 0
  w <- drop((1 - W) * exp(Ze_Beta[[1]] %*% theta))
  a <- matrix(rep(0, p_Z^2), nrow = p_Z, ncol = p_Z)
  for (j in 1:N) {
    a <- a + (Ze_Beta[[1]][j,]) %*% t(Ze_Beta[[1]][j,]) * w[j]
  }
  w2 <- a * sum(w) - t(w %*% Ze_Beta[[1]]) %*% (w %*% Ze_Beta[[1]])
  f <- f + drop(w2 / (sum(w))^2)
  return(f)
}

# ======================================== #
#              Method of CEB               #
# ======================================== #

# CEB / EBME: our proposed method----
objective.ebme.r_Beta <- function(theta) {
  f <- drop(log((1 - W) %*% exp(Ze_Beta[[1]] %*% theta)) -
              Ze1_bar_Beta[[1]] %*% theta -
              t(theta) %*% Sigma %*% theta / 2)
  return(f)
}

gradient.ebme.r_Beta <- function(theta) {
  w <- drop((1 - W) * exp(Ze_Beta[[1]] %*% theta))
  f <- drop(w %*% Ze_Beta[[1]] / sum(w) - Ze1_bar_Beta[[1]]) - drop(Sigma %*% theta)
  return(f)
}

# ======================================== #
#              Method of HL                #
# ======================================== #

# Chengcheng Hu & D. Y Lin (2004): replicate method 1----
gradient.hcc.r_Beta <- function(theta) {
  eta0 <- 0
  eta1 <- 0
  for (u in 2:r) {
    for (s in 1:(u - 1)) {
      eta0 <- eta0 + sum(exp(drop((theta %*% t(Ze_Beta[[u]] - Ze_Beta[[s]])))) +
                           exp(drop((theta %*% t(Ze_Beta[[s]] - Ze_Beta[[u]])))))
      eta1 <- eta1 +
        (drop(t(Ze_Beta[[u]] - Ze_Beta[[s]]) %*% (exp(drop(theta %*% t(Ze_Beta[[u]] - Ze_Beta[[s]])))))) +
        (drop(t(Ze_Beta[[s]] - Ze_Beta[[u]]) %*% (exp(drop(theta %*% t(Ze_Beta[[s]] - Ze_Beta[[u]]))))))
    }
  }
  eta0 <- sqrt(eta0 / (r * (r - 1) * N))
  eta1 <- eta1 / (2 * N * r * (r - 1) * eta0)
  f0 <- 0
  f1 <- 0
  for (k in 1:r) {
    f0 <- f0 + sum(drop((1 - W) * exp(Ze_Beta[[k]] %*% theta)))
    f1 <- f1 + drop((t(Ze_Beta[[k]]) - eta1 / eta0) %*% drop((1 - W) * exp(Ze_Beta[[k]] %*% theta)))
  }
  f <- drop(f1 / f0 - Ze1_bar_al_Beta)
  # f <- sum(f^2)
  return(f)
}


# ======================================== #
#              Method of HW                #
# ======================================== #
# - Only the estimating equation for mismeasured covariates are constructed with the method of Huang & Wang.
# - The estimating equation for precisely measuered covariates are constructed with the method of Hu & Lin.

# Yijian Huang & C. Y. Wang (2000): replicate method 2----
gradient.hyj.r_Beta <- function(theta) {

  index_Permutation <- permutations(n = r, r = 2, v = 1:r)
  N_Permutation <- nrow(index_Permutation)
  # RepDataset_Z <- list(repset1 = Ze_Beta[[1]], repset2 = Ze_Beta[[2]], repset3 = Ze_Beta[[3]])
  RepDataset_Z <- list(repset1 = Ze_Beta[[1]], repset2 = Ze_Beta[[2]])

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
      Denominator <- (1 - W) %*% exp(Ze_Beta[[k]] %*% theta)
    }else {
      Denominator <- Denominator + (1 - W) %*% exp(Ze_Beta[[k]] %*% theta)
    }
  }
  Denominator <- Denominator / r
  # 注意：后半部分Z1bar【不能】直接用colMeans！m1要除的是n1，不是N！
  HYJ_val <- Numerator / drop(Denominator) - Ze1_bar_al_Beta
  return(HYJ_val)
}
