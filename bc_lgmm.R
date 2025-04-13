

Optbw.bclgmm <- function(Xmat, Y, Zmat, gmmwl, kernel_type = "epanechnikov") {
  bw_list <- seq(0.4, 0.81, by = 0.03)
  T <- nrow(Y)
  ssr <- sapply(bw_list, function(h) {
    est <- tryCatch(lgmm2_loo(X = Xmat, Y =Y, Z = Zmat, bw = T^h, 
                              gmmwl = gmmwl, 
                              kernel_type = "epanechnikov"), 
                    error = function(e) return(NULL))
    if (!is.null(est)) sum(unlist(est$res)^2) else NA
  })
  valid <- na.omit(data.frame(bw = bw_list, SSR = ssr))
  valid <- valid[order(valid$SSR), ]
  return(valid[which.min(valid[, 2]), 1])
}


lgmm2_loo <- function(X, Y, Z, bw, gmmwl,
                      kernel_type = "epanechnikov") {
  nobs <- nrow(Y)
  c <- matrix(1, nobs, 1)

  k <- ncol(Z)
  
  bw0 <- floor(bw)
  Xmat_mir <- mirror_data(X, bw0)
  Zmat_mir <- mirror_data(Z, bw0)
  Y_mir <- mirror_data(Y, bw0)
  
  p <- ncol(X)
  Lt <- vector("list", nobs)
  Ht <- vector("list", nobs)
  beta_mat <- matrix(NA, nrow = p, ncol = nobs)
  res <- vector("list", nobs)
  
  
  ## === 1) Local Regression Stage ===
  for (new_t in 1:nobs) {
    t_new <- new_t + bw0
    Z_loo <- Zmat_mir[-t_new, , drop = FALSE]
    XS_loo <- Xmat_mir[-t_new, , drop = FALSE]
    Y_loo <- Y_mir[-t_new, , drop = FALSE]
    

    K_t <- kernelt(t_new, (nobs + 2*bw0), kernel_type, bandwidth = bw)
    K_t <- K_t[-t_new,-t_new]
    Gt <- t(Z_loo) %*% K_t %*% XS_loo
    ht <- t(Z_loo) %*% K_t %*% Y_loo
    beta_mat[, new_t] <- solve(t(Gt) %*% solve(gmmwl[[new_t]]) %*% Gt) %*%
      t(Gt) %*% solve(gmmwl[[new_t]]) %*% ht
    
    
    res[[new_t]] <- Y_mir[new_t] - as.numeric(Xmat_mir[new_t, ] %*% t(beta_mat[, new_t]))
    
  }
  
  list(res = res)
  }
  





lgmm2 <- function(X, Y, Z, W, bw = NULL, gmmw = NULL, 
                  withZ = 1, lag = 0, kernel_type = "epanechnikov") {
  nobs <- nrow(Y)
  c <- matrix(1, nobs, 1)
  if (withZ == 1) {
    Xmat <- as.matrix(data.frame(c, X, Z))
    Zmat <- as.matrix(data.frame(c, W, Z))
  } else {
    Xmat <- as.matrix(data.frame(c, X))
    Zmat <- as.matrix(data.frame(c, W))
  }
  
  k <- ncol(Zmat)
  if (is.null(gmmw)) {
    cat("gmmw is NULL — using identity matrices.\n")
    gmmwl <- replicate(nobs, diag(k), simplify = FALSE)
  } else {
    gmmwl <- gmmw
  }

    
  if (is.null(bw) || is.na(bw)) {
    opt_bw <- Optbw.bclgmm(Xmat = Xmat, Y = Y, Zmat=Zmat, gmmwl = gmmwl, kernel_type = "epanechnikov")
    bw <- nobs^opt_bw
    print(paste("Optimal bandwidth selected:", opt_bw))
  }

  bw0 <- floor(bw)
  Xmat_mir <- mirror_data(Xmat, bw0)
  Zmat_mir <- mirror_data(Zmat, bw0)
  Y_mir    <- mirror_data(Y, bw0)
  
  p <- ncol(Xmat)
  Lt <- vector("list", nobs)
  Ht <- vector("list", nobs)
  beta_mat <- matrix(NA, nrow = p, ncol = nobs)
  
  ## === 1) Local Regression Stage ===
  for (new_t in 1:nobs) {
    t_new <- new_t + bw0
    K_t <- kernelt(t_new, (nobs + 2*bw0), kernel_type, bandwidth = bw)
    print(new_t)
    Gt <- t(Zmat_mir) %*% K_t %*% Xmat_mir
    ht <- t(Zmat_mir) %*% K_t %*% Y_mir
    beta_mat[, new_t] <- solve(t(Gt) %*% solve(gmmwl[[new_t]]) %*% Gt) %*%
      t(Gt) %*% solve(gmmwl[[new_t]]) %*% ht
    
    Lt[[new_t]] <- - t(Zmat_mir) %*% K_t %*% Xmat_mir
    Ht[[new_t]] <- solve(t(Lt[[new_t]]) %*% solve(gmmwl[[new_t]]) %*% Lt[[new_t]]) %*%
      t(Lt[[new_t]]) %*% solve(gmmwl[[new_t]])
  }
  
  beta_mat_mirr <- t(mirror_data(t(beta_mat), bw0))
  
  ## === 2) Local HAC Variance Part ===
  M2l <- vector("list", nobs)
  for (new_t in 1:nobs) {
    t_new <- new_t + bw0
    K_t <- kernelt(t_new, (nobs + 2*bw0), kernel_type, bandwidth = bw0)
    
    s_start <- t_new - bw0
    s_end   <- t_new + bw0
    M2 <- matrix(0, nrow = k, ncol = k)
    for (s in s_start:s_end) {
      k_st <- diag(K_t)[s]
      res_s <- Y_mir[s] - as.numeric(Xmat_mir[s, ] %*% beta_mat_mirr[, s])
      m_s <- as.matrix(Zmat_mir[s, ]) * res_s
      M2 <- M2 + (k_st^2) * (m_s %*% t(m_s))
    }
    M2l[[new_t]] <- M2
  }
  
  M_crossl <- vector("list", nobs)
  M_suml   <- vector("list", nobs)
  Q_l      <- vector("list", nobs)
  for (new_t in 1:nobs) {
    t_new <- new_t + bw0
    K_t <- kernelt(t_new, (nobs + 2*bw0), kernel_type, bandwidth = bw0)
    
    M_cross <- matrix(0, nrow = k, ncol = k)
    for (j in (1:(2 * bw0))) {
      w_j <- bartlett_kernel(j, L = lag)
      s_start <- t_new - bw0
      s_end   <- t_new + bw0 - j
      if (s_end >= s_start & (lag != 0)) {
        for (s in s_start:s_end) {
          m_s  <- as.matrix(Zmat_mir[s, ]) *
            (Y_mir[s] - as.numeric(Xmat_mir[s, ] %*% beta_mat_mirr[, s]))
          m_sj <- as.matrix(Zmat_mir[s + j, ]) *
            (Y_mir[s + j] - as.numeric(Xmat_mir[s + j, ] %*% beta_mat_mirr[, s + j]))
          M_cross <- M_cross +
            w_j * diag(K_t)[s] * diag(K_t)[s+j] * (m_s %*% t(m_sj))
        }
      }
    }
    M_crossl[[new_t]] <- (2 * M_cross) 
    M_suml[[new_t]]   <- M2l[[new_t]] + M_crossl[[new_t]]
    Q_l[[new_t]]      <- Ht[[new_t]] %*% M_suml[[new_t]] %*% t(Ht[[new_t]])
  }
  
  list(beta = beta_mat, var = Q_l, opt_gmmw = M_suml, bw = bw0)
}





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

g <- stage1$opt_gmmw
bw2 <- stage1$bw
stage2 <- lgmm2(
  X = X, Y = Y, Z = Z, W = W, bw = bw2,
  gmmw = g, withZ = 0,
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

library(ggplot2)
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







