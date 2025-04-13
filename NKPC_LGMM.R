library(zoo)
library(ggplot2)
setwd("E:/Spring 2025/Research/Bertille/Code/")



data <- read.csv("App_Data/data_main_complete.csv")  

data$obs <- as.Date(as.yearmon(data$obs, "%Y M%m"))


data <- na.omit(data)

Y <- as.matrix(data$d_I)
X <- as.matrix(data[, c("d_UNEMP")])
W <- as.matrix(data[, c( "d_UNEMP_1", "d_UNEMP_2", "d_UNEMP_3", "d_UNEMP_4")])
Z <- as.matrix(data[,c("d_I_1")])



stage1 <- lgmm2(X, Y, Z, W,
      bw = NULL, gmmw = NULL, 
      withZ = 1, lag = 0,
      kernel_type = "epanechnikov")


g <- stage1$opt_gmmw
bw2 <- stage1$bw
stage2 <- lgmm2(X, Y, Z, W,
      bw = bw2, gmmw = g, 
      withZ = 1, lag = 0,
      kernel_type = "epanechnikov")


saveRDS(stage1, file = "stage1_bclgmm.rds")
saveRDS(stage2, file = "stage2_bclgmm.rds")

beta_X <- stage2$beta[2, ]
se_X   <- sapply(stage2$var, function(v) sqrt(v[2, 2]))
lower  <- beta_X - 1.96 * se_X
upper  <- beta_X + 1.96 * se_X

plot_df <- data.frame(
  Time       = data$obs,
  theta_hat  = beta_X,
  lower      = lower,
  upper      = upper
)


ggplot(plot_df, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = theta_hat), color = "blue", size = 0.6) +
  labs(
    title = "LGMM Estimation of Inflation Change Coefficient",
    y = expression(Coefficient),
    x = expression(Time)
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.4))

