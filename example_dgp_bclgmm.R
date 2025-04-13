library(ggplot2)

source("lgmm_paper_functions.R")



# ========== Example Usage ==========
set.seed(123)
T <- 100

data <- dgp_e2(T)  

X <- as.matrix(data$X_t)
Y <- as.matrix(data$Y_t)
W <- as.matrix(data[, c("Z_t1", "Z_t2", "Z_t3", "Z_t4")])
Z <- matrix(0, nrow = T, ncol = 0)


stage1 <- lgmm2(
  X = X, Y = Y, Z = Z, W = W, bw = NULL,
  gmmw = NULL, withZ = 0,
  lag = 0, kernel_type = "epanechnikov"
)

gmw <- stage1$opt_gmmw
bw2 <- stage1$bw

stage2 <- lgmm2(
  X = X, Y = Y, Z = Z, W = W, bw = bw2,
  gmmw = gmw, withZ = 0,
  lag = 0, kernel_type = "epanechnikov"
)

beta_X <- stage2$beta[2, ]
se_X   <- sapply(stage2$var, function(v) sqrt(v[2, 2]))
lower  <- beta_X - 1.96 * se_X
upper  <- beta_X + 1.96 * se_X

plot_df <- data.frame(
  Time       = 1:T,
  tau        = data$tau,
  theta_true = data$theta_2,
  theta_hat  = beta_X,
  lower      = lower,
  upper      = upper
)


ggplot(plot_df, aes(x = tau)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = theta_hat), color = "blue", size = 0.6) +
  geom_line(aes(y = theta_true), color = "red", linetype = "dashed", size = 0.6) +
  labs(
    title = "LGMM Estimation of Time-Varying Coefficient (theta_2)",
    y = expression(theta[2](tau)),
    x = expression(tau)
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
