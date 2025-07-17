#-----------------------------------------------
# Simulation: Impact of Multicollinearity & Measurement Error on Estimator Performance
# Outputs: Excel results, convergence vs CI plots (plain & labeled), combined bias/MSE facet plots (plain & labeled)
# 相比v2，增加结果表格输出
#-----------------------------------------------

# 0. Setup
tmp <- ls()
rm(list = tmp)
options(warn = -1)

library(MASS)      # mvrnorm
library(gtools)    # permutations()
library(car)       # vif()
library(ggplot2)   # plotting
library(patchwork) # plot composition
library(rootSolve) # multiroot
library(openxlsx)  # Excel export
library(dplyr)
library(tidyr)

# load objective/gradient/hessian
source("./CEB/R/09SuppSim_Multicollinearity/999.ObjectiveGradient_NormalME.R")

# 1. Parameters
set.seed(102)
N <- 1000
N_Sim <- 1000
r <- 2
Mean_X <- c(4, 2)
Mean_U <- c(3, 1)
p_X <- length(Mean_X)
p_U <- length(Mean_U)
p_Z <- p_X + p_U
# vec_rhoZ <- c(0, seq(0.1, 0.5, by = 0.2))
vec_rhoZ <- seq(0.0, 0.9, by = 0.3)
vec_eps <- seq(0.1, 0.5, by = 0.2)
True.theta <- c(3.5, -1, 0.5, -0.25, -0.1)
Meancoef <- c(210, 27.4, 13.7, 13.7, 13.7)
Sigma_Y <- 4
rho_Y <- 0
ATE <- 10

# 2. Results container
df_res <- expand.grid(rhoZ = vec_rhoZ, eps = vec_eps,
                      method = c("EB", "CEB", "BCEB", "CEB-HL", "CEB-HW"),
                      stringsAsFactors = FALSE)
df_res <- mutate(df_res, bias = NA, mse = NA, CI = NA, conv_rate = NA)

# styles
shapes <- c(EB = 16, CEB = 17, BCEB = 15, `CEB-HL` = 3, `CEB-HW` = 0)
colors <- c(EB = "#FC8D62", CEB = "#8DA62C", BCEB = "#66C2A5",
            `CEB-HL` = "#8DA0CB", `CEB-HW` = "#E78AC3")
# linetypes <- c(EB=1, CEB=2, BCEB=3, `CEB-HL`=4, `CEB-HW`=5)

# 3. Simulation with progress bar
total <- length(vec_rhoZ) * length(vec_eps) * N_Sim
pb <- txtProgressBar(min = 0, max = total, style = 3)
count <- 0
for (rhoZ in vec_rhoZ) {
  Sigma_Z <- matrix(rhoZ, nrow = p_Z, ncol = p_Z)
  diag(Sigma_Z) <- 1
  for (eps in vec_eps) {
    Sigma_e <- diag(rep(eps, p_X))
    Sigma <- diag(c(diag(Sigma_e), rep(0, p_U)))
    # convergence flags
    iter_flags <- matrix(FALSE, nrow = N_Sim, ncol = 5)
    colnames(iter_flags) <- c("EB", "CEB", "BCEB", "CEB-HL", "CEB-HW")
    ATT_mat <- matrix(NA, nrow = N_Sim, ncol = 5, dimnames = list(NULL, colnames(iter_flags)))
    CI_vec <- numeric(N_Sim)
    for (i in 1:N_Sim) {
      # data
      Z_true <- mvrnorm(N, mu = c(Mean_X, Mean_U), Sigma = Sigma_Z)
      X_reps <- lapply(1:r, function(k) Z_true[, 1:p_X] + mvrnorm(N, rep(0, p_X), Sigma_e))
      Ze_list <- lapply(X_reps, function(Xr) {
        mat <- cbind(Xr, Z_true[, (p_X + 1):p_Z])
        colnames(mat) <- paste0("V", 1:p_Z)
        mat
      })
      assign("Ze_Normal", Ze_list, envir = .GlobalEnv)
      W <- rbinom(N, 1, plogis(cbind(1, Z_true) %*% True.theta))
      # assign("W", W, envir = .GlobalEnv)
      # assign("Sigma", Sigma, envir = .GlobalEnv)
      # assign("r", r, envir = .GlobalEnv)
      # assign("p_Z", p_Z, envir = .GlobalEnv)
      # assign("N", N, envir = .GlobalEnv)
      idx_t <- which(W == 1)
      idx_c <- which(W == 0)
      Ze0 <- Ze_list[[1]][idx_c, , drop = FALSE]
      Ze1_bars <- lapply(Ze_list, function(ZL)colMeans(ZL[idx_t, , drop = FALSE]))
      assign("Ze1_bar_Normal", Ze1_bars, envir = .GlobalEnv)
      assign("Ze1_bar_al_Normal", Reduce(`+`, Ze1_bars) / r, envir = .GlobalEnv)
      # CI
      df_obs <- as.data.frame(cbind(Reduce(`+`, X_reps) / r, Z_true[, (p_X + 1):p_Z]))
      eigs <- eigen(t(scale(df_obs, TRUE, FALSE)) %*% scale(df_obs, TRUE, FALSE))$values
      CI_vec[i] <- sqrt(max(eigs) / min(eigs))
      # outcomes
      # treatment effect
      Mean_i_0 <- cbind(1, Z_true) %*% Meancoef
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
      Y1 <- Y[idx_t]
      Y0 <- Y[idx_c]

      # EB: naive entropy balancing
      res1 <- optim(par = rep(0, p_Z),
                    fn = objective.zqy.r_Normal,
                    gr = gradient.zqy.r_Normal,
                    method = "BFGS",
                    control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200))
      th1 <- res1$par
      iter_flags[i, 'EB'] <- (res1$convergence == 0)
      w1 <- exp(Ze0 %*% th1)
      w1 <- w1 / sum(w1)
      ATT_mat[i, "EB"] <- mean(Y1) - sum(w1 * Y0)

      # CEB: entropy balancing with measurement-error penalty
      output.CEB.solveObj_Normal <- tryCatch(optim(par = rep(0, p_Z),
                                                   fn = objective.ebme.r_Normal,
                                                   gr = gradient.ebme.r_Normal,
                                                   method = "BFGS",
                                                   control = list(trace = 0, reltol = .Machine$double.eps, maxit = 200)),
                                             error = function(e) {
                                               print("Fail to solve objective of CEB, We directly use naive EB's results for subsequent comparison.")
                                               # 此时无法求解CEB，直接将solve CEB Objective的结果赋为FALSE
                                               return(FALSE)
                                             })
      # 如果optim无法求解，output.CEB.solveObj_Normal会是一个logical变量；否则，如果optim能求解，output.CEB.solveObj_Normal是一个list变量
      if (class(output.CEB.solveObj_Normal) != "list" |
        is.na(output.CEB.solveObj_Normal$value) |
        is.nan(output.CEB.solveObj_Normal$value)) {
        # 此时无法求解CEB，将naive EB的结果赋给CEB
        output.CEB.solveObj_Normal <- res1
      }
      theta_CEB_solveObj <- output.CEB.solveObj_Normal$par
      ObjValue_CEB_solveObj <- output.CEB.solveObj_Normal$value

      # theta_CEB - solve the gradient funtion
      output.CEB.solveGradient_Normal <- tryCatch(multiroot(f = gradient.ebme.r_Normal,
                                                            start = rep(0, p_Z),
                                                            rtol = .Machine$double.eps, maxiter = 200),
                                                  error = function(e) {
                                                    print("Fail to solve gradient of CEB, We directly use naive EB's results for subsequent comparison.")
                                                    # 此时无法求解CEB，直接将solve CEB Objective的结果赋为FALSE
                                                    return(FALSE)
                                                  })
      # 如果multiroot无法求解，output.CEB.solveGradient会是一个logical变量；否则，如果multiroot能求解，output.CEB.solveGradient是一个list变量
      if (class(output.CEB.solveGradient_Normal) != "list") {
        # 此时求解CEB gradient 不收敛，直接赋值9999
        ObjValue_CEB_solveGradient <- 9999
      }else {
        theta_CEB_solveGradient <- output.CEB.solveGradient_Normal$root
        ObjValue_CEB_solveGradient <- objective.ebme.r_Normal(theta = theta_CEB_solveGradient)
        if (is.nan(ObjValue_CEB_solveGradient) | is.na(ObjValue_CEB_solveGradient)) { ObjValue_CEB_solveGradient <- 9999 }
      }

      if (ObjValue_CEB_solveObj <= ObjValue_CEB_solveGradient) {
        th2 <- theta_CEB_solveObj
        iter_flags[i, 'CEB'] <- (max(abs(gradient.ebme.r_Normal(theta_CEB_solveObj))) < 1e-04)
        ## double check
        if (is.na(iter_flags[i, 'CEB']) | is.nan(iter_flags[i, 'CEB'])) {
          iter_flags[i, 'CEB'] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }else {
        th2 <- theta_CEB_solveGradient
        iter_flags[i, 'CEB'] <- (max(abs(output.CEB.solveGradient_Normal$f.root)) < sqrt(.Machine$double.eps)) & (output.CEB.solveGradient_Normal$estim.precis != "NaN")
        ## double check
        if (is.na(iter_flags[i, 'CEB']) | is.nan(iter_flags[i, 'CEB'])) {
          iter_flags[i, 'CEB'] <- FALSE
          print("Solution to CEB fail to converge!")
        }
      }
      w2 <- exp(Ze0 %*% th2)
      w2 <- w2 / sum(w2)
      ATT_mat[i, "CEB"] <- mean(Y1) - sum(w2 * Y0)

      # BCEB: bias-corrected EB
      H <- hessian.zqy.r_Normal(th1)
      th3 <- tryCatch(solve(H - Sigma, H %*% th1), error = function(e) rep(NaN, p_Z))
      iter_flags[i, 'BCEB'] <- all(is.finite(th3))
      w3 <- exp(Ze0 %*% th3)
      w3 <- w3 / sum(w3)
      ATT_mat[i, "BCEB"] <- mean(Y1) - sum(w3 * Y0)

      # CEB-HL: Huang & Lin (HCC)
      sol4 <- tryCatch(multiroot(gradient.hcc.r_Normal, rep(0, p_Z)), error = function(e) NULL)
      if (is.null(sol4)) {
        iter_flags[i, 'CEB-HL'] <- FALSE
        th4 <- rep(NA, p_Z)
      } else {
        iter_flags[i, 'CEB-HL'] <- TRUE
        th4 <- sol4$root
      }
      w4 <- exp(Ze0 %*% th4)
      w4 <- w4 / sum(w4)
      ATT_mat[i, "CEB-HL"] <- mean(Y1) - sum(w4 * Y0)
      ### double check
      if (is.na(iter_flags[i, 'CEB-HL'])) {
        iter_flags[i, 'CEB-HL'] <- FALSE
        print("NA clear!")
      }
      ### double check
      if (is.na(ATT_mat[i, "CEB-HL"])) {
        iter_flags[i, 'CEB-HL'] <- FALSE
        print("ATT of HL - NA clear!")
      }

      # CEB-HW: Huang & Wang (HYJ)
      sol5 <- tryCatch(multiroot(gradient.hyj.r_Normal, rep(0, p_Z)), error = function(e) NULL)
      if (is.null(sol5)) {
        iter_flags[i, 'CEB-HW'] <- FALSE
        th5 <- rep(NA, p_Z)
      } else {
        iter_flags[i, 'CEB-HW'] <- TRUE
        th5 <- sol5$root
      }
      iter_flags[i, 'CEB-HW'] <- (max(gradient.hyj.r_Normal(sol5$root)) < sqrt(.Machine$double.eps)) & (sol5$estim.precis != "NaN")
      w5 <- exp(Ze0 %*% th5)
      w5 <- w5 / sum(w5)
      ATT_mat[i, "CEB-HW"] <- mean(Y1) - sum(w5 * Y0)
      ### double check
      if (is.na(iter_flags[i, 'CEB-HW'])) {
        iter_flags[i, 'CEB-HW'] <- FALSE
        print("NA clear!")
      }
      ### double check
      if (is.na(ATT_mat[i, "CEB-HW"])) {
        iter_flags[i, 'CEB-HW'] <- FALSE
        print("ATT of HW - NA clear!")
      }

      # aggregate subresults
      sel <- df_res$rhoZ == rhoZ & df_res$eps == eps
      df_res$CI[sel] <- mean(CI_vec[1:i])
      for (meth in df_res$method[sel]) {
        idxc <- iter_flags[1:i, meth]
        vals <- ATT_mat[idxc, meth]
        df_res$bias[sel & df_res$method == meth] <- mean(vals - ATE, na.rm = TRUE)
        df_res$mse[sel & df_res$method == meth] <- mean((vals - ATE)^2, na.rm = TRUE)
        df_res$conv_rate[sel & df_res$method == meth] <- mean(idxc)
      }
      count <- count + 1
      setTxtProgressBar(pb, count)
    }
  }
}
close(pb)

# 4. Export to Excel
write.xlsx(df_res, file = paste0("BiasMSEofATT_Rc_CI_N", N / 1000, "k_Sim", N_Sim / 1000, "k_rhoZ.xlsx"), rowNames = FALSE)

# 5a. Convergence vs rhoZ: plain & labeled
# 1) 构造 conv_df
conv_df <- df_res %>%
  group_by(
    eps,
    rhoZ,
    method
  ) %>%
  summarize(
    # CI = mean(CI, na.rm = TRUE),
    conv_rate = mean(conv_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    eps_lab = paste0("Var(epsilon)==", eps),
    method = factor(method, levels = c("EB", "CEB", "BCEB", "CEB-HL", "CEB-HW"))
  )
# 2) 画图
p_conv_facet <- ggplot(conv_df,
                       aes(x = rhoZ,
                           y = conv_rate,
                           color = method,
                           shape = method,
                           group = method)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~eps_lab,
             scales = "free_x",         # 每个子图横轴自适应
             nrow = 1,                # 一行排列
             labeller = label_parsed      # 解析 Var(ε)==0.1
  ) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  # scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_x_continuous(
    breaks = vec_rhoZ
    # , # 你自己的刻度
    # labels = function(x) parse(text = paste0("rho[Z]==", x))
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(5)) +
  labs(
    x = expression(rho[Z]),
    y = "Convergence Rate",
    color = "Method",
    shape = "Method") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# # labeled y-values
# p_conv_label <- p_conv_plain + geom_text(aes(label=round(conv,2)),vjust=-0.5,size=3)
# 3) 保存
ggsave("convergence_vs_rhoZ_by_eps_rhoZ.pdf",
       p_conv_facet,
       width = 10,
       height = 4)
# ggsave("convergence_vs_CI_labeled.pdf",p_conv_label,width=6,height=5)

# 5b. Bias/MSE facet: plain & labeled
agg2 <- df_res %>%
  group_by(rhoZ, eps, method) %>%
  summarize(bias = mean(bias), mse = mean(mse),
            # CI = mean(CI),
            .groups = 'drop')
df_long <- agg2 %>%
  pivot_longer(c(bias, mse), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("bias", "mse"), labels = c("Bias", "MSE")),
         eps_lab = paste0("Var(epsilon)==", eps),
         method = factor(method, levels = names(shapes)))
common_aes2 <- aes(
  # x = CI,
  x = rhoZ,
  y = Value, color = method, shape = method, group = method)
# plain
p_comb_plain <- ggplot(df_long, common_aes2) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(values = colors, breaks = names(shapes)) +
  scale_shape_manual(values = shapes, breaks = names(shapes)) +
  # scale_linetype_manual(values=linetypes,breaks=names(shapes)) +
  scale_x_continuous(
    breaks = vec_rhoZ
    # , # 你自己的刻度
    # labels = function(x) parse(text = paste0("rho[Z]==", x))
  ) +
  facet_grid(Metric ~ eps_lab, scales = "free_y", switch = "y", labeller = label_parsed) +
  labs(
    # x = "Condition Index (CI)",
    x = expression(rho[Z]),
    y = NULL, color = "Method", shape = "Method") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey80", color = NA), strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90"), legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
# # labeled values
# p_comb_label <- p_comb_plain + geom_text(aes(label=round(Value,2)),vjust=-0.5,size=2)
ggsave(paste0("BiasMSEofATTvsrhoZ_plain_N", N / 1000, "k_Sim", N_Sim / 1000, "k_rhoZ.pdf"),
       p_comb_plain,
       width = 10, height = 6)
# ggsave("combined_facet_labeled.pdf",p_comb_label,width=10,height=6)

save.image("998.CI_BiasMSEofATT_RateConvergence_VecrhoZv3.Rdata")