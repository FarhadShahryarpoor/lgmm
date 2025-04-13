dgp_e2 <- function(T) {
  rho <- 0.8 
  tau <- (1:T) / T 
  theta_1_tau <- 0.2 * exp(-0.7 + 3.5 * tau)
  theta_2_tau <- 2 * tau + exp(-16 * (tau - 0.5)^2) - 1
  
  u_t <- mvrnorm(T, mu = rep(0, 4), Sigma = diag(4))
  u_t1 <- u_t[, 1]
  u_t2 <- u_t[, 2]
  u_t3 <- u_t[, 3]
  u_t4 <- u_t[, 4]
  
  Z_t1 <- numeric(T)
  Z_t2 <- numeric(T)
  Z_t3 <- numeric(T)
  Z_t4 <- numeric(T)
  
  Z_t1[1] <- 0
  Z_t2[1] <- 0
  Z_t3[1] <- 0
  Z_t4[1] <- 0
  
  for (t in 2:T) {
    Z_t1[t] <- 0.9 * sin(4 * pi * tau[t]) * Z_t1[t - 1] + u_t1[t]
    Z_t2[t] <- 0.5 + 0.8 * sin(4 * pi * tau[t]) * Z_t2[t - 1] + u_t2[t]
    Z_t3[t] <- 0.9 * sin(4 * pi * tau[t]) * Z_t3[t - 1] + u_t3[t]
    Z_t4[t] <- 0.9 * sin(4 * pi * tau[t]) * Z_t4[t - 1] + u_t4[t]
  }
  
  Z_t <- cbind(Z_t1, Z_t2, Z_t3, Z_t4)
  
  errors <- mvrnorm(T, mu = c(0, 0), Sigma = matrix(c(1, rho, rho, 1), 2, 2))
  e_t <- errors[, 1]
  nu_t <- errors[, 2]
  
  X_t <- 0.8 * Z_t[, 1] + 1.3 * Z_t[, 2] + Z_t[, 3] + Z_t[, 4] + nu_t
  
  Y_t <- theta_1_tau + theta_2_tau * X_t + e_t
  
  
  cons <- rep(1,T)
  data.frame(tau = tau, X_t = X_t, Y_t = Y_t, cons, Z_t1 = Z_t1, Z_t2 = Z_t2, Z_t3 = Z_t3, Z_t4 = Z_t4, theta_1 = theta_1_tau, theta_2 = theta_2_tau)
}




bartlett_kernel <- function(j, L) {
  j_abs <- abs(j)
  if (j_abs <= L) {
    return(1 - (j_abs / (L + 1)))
  } else {
    return(0)
  }
}


kernelt <- function(t, nobs, kernel_type = "gaussian", bandwidth = 5) {
  if (kernel_type == "gaussian") {
    kernel_function <- function(s, t) {
      exp(- (s - t)^2 / (2 * bandwidth^2))
    }
  } else if (kernel_type == "uniform") {
    kernel_function <- function(s, t) {
      ifelse(abs(s - t) <= bandwidth, 1, 0)
    }
  } else if (kernel_type == "triangular") {
    kernel_function <- function(s, t) {
      w <- 1 - abs(s - t) / bandwidth
      ifelse(abs(s - t) <= bandwidth, w, 0)
    }
  } else if (kernel_type == "epanechnikov") {
    kernel_function <- function(s, t) {
      u <- (s - t) / bandwidth
      ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
    }
  } else {
    stop("Kernel type not supported. Choose from 'gaussian', 'uniform', 'triangular', or 'epanechnikov'.")
  }
  
  kernel_weights <- sapply(1:nobs, function(s) kernel_function(s, t))
  K_t <- diag(kernel_weights)
  return(K_t)
}


Optbw.lgmm <- function(Xmat, Y, Zmat, gmmwl, kernel_type = "epanechnikov") {
  bw_list <- seq(0.4, 0.81, by = 0.03)
  T <- nrow(Y)
  ssr <- sapply(bw_list, function(h) {
    est <- tryCatch(lgmm_loo(X = Xmat, Y =Y, Z = Zmat, bw = T^h, 
                              gmmwl = gmmwl, 
                              kernel_type = "epanechnikov"), 
                    error = function(e) return(NULL))
    if (!is.null(est)) sum(unlist(est$res)^2) else NA
  })
  valid <- na.omit(data.frame(bw = bw_list, SSR = ssr))
  valid <- valid[order(valid$SSR), ]
  return(valid[which.min(valid[, 2]), 1])
}


lgmm_loo <- function(X, Y, Z, bw, gmmwl,
                     kernel_type = "epanechnikov"){
  nobs <- nrow(Y)
  c <- matrix(1, nobs, 1)

  k <- ncol(Z)
  
  p <- ncol(X)
  Lt <- vector("list", nobs)
  Ht <- vector("list", nobs)
  beta_mat <- matrix(NA, nrow = p, ncol = nobs)
  res <- vector("list", nobs)
  
  ## =========== 1) Local Regression Stage =============
  for (t in 1:nobs) {
    
    Z_loo <- Z[-t, , drop = FALSE]
    XS_loo <- X[-t, , drop = FALSE]
    Y_loo <- Y[-t, , drop = FALSE]
    
    
    K_t <- kernelt(t, nobs, kernel_type = kernel_type, bandwidth = bw)
    K_t <- K_t[-t,-t]
    Gt <- t(Z_loo) %*% K_t %*% XS_loo
    ht <- t(Z_loo) %*% K_t %*% Y_loo
    
    beta_mat[, t] <-
      solve(t(Gt) %*% solve(gmmwl[[t]]) %*% Gt) %*%
      t(Gt) %*% solve(gmmwl[[t]]) %*% ht
    
    res[[t]] <- Y[t] - as.numeric(X[t, ] %*% t(beta_mat[, t]))
  }
  list(
    res = res
  )
}




lgmm <- function(X, Y, Z, W, bw = NULL, gmmw = NULL,
                 withZ = 1, lag = 0, kernel_type = "epanechnikov")
{
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
    opt_bw <- Optbw.lgmm(Xmat = Xmat, Y = Y, Zmat=Zmat, gmmwl = gmmwl, kernel_type = "epanechnikov")
    bw <- floor(nobs^opt_bw)
    print(paste("Optimal bandwidth selected:", opt_bw))
  }
  
  
  p <- ncol(Xmat)
  Lt <- vector("list", nobs)
  Ht <- vector("list", nobs)
  beta_mat <- matrix(NA, nrow = p, ncol = nobs)
  
  ## =========== 1) Local Regression Stage =============
  for (t in 1:nobs) {
    K_t <- kernelt(t, nobs, kernel_type = kernel_type, bandwidth = bw)
    
    Gt <- t(Zmat) %*% K_t %*% Xmat
    ht <- t(Zmat) %*% K_t %*% Y
    
    beta_mat[, t] <-
      solve(t(Gt) %*% solve(gmmwl[[t]]) %*% Gt) %*%
      t(Gt) %*% solve(gmmwl[[t]]) %*% ht
    
    Lt[[t]] <- -1 * t(Zmat) %*% K_t %*% Xmat
    Ht[[t]] <-
      solve(t(Lt[[t]]) %*% solve(gmmwl[[t]]) %*% Lt[[t]]) %*%
      t(Lt[[t]]) %*% solve(gmmwl[[t]]) #H_t matrix for the breads of VC matrix
  }
  
  ## =========== 2) Local HAC Variance Part =============
  M2l <- vector("list", nobs)
  for (t in 1:nobs) {
    K_t <- kernelt(t, nobs, kernel_type = kernel_type, bandwidth = bw)
    
    s_start <- max(1, t - floor(bw))
    s_end   <- min(nobs, t + floor(bw))
    
    M2 <- matrix(0, nrow = ncol(Zmat), ncol = ncol(Zmat))
    
    for (s in s_start:s_end) {
      k_st <- diag(K_t)[s]
      w2 <- k_st^2
      res_s <- Y[s] - as.numeric(Xmat[s, ] %*% beta_mat[, s])
      m_s   <- as.matrix(Zmat[s, ]) * res_s
      M2    <- M2 + w2 * (m_s %*% t(m_s))
    }
    M2l[[t]] <- M2 
  }
  
  M_crossl <- vector("list", nobs)
  Q_l <- vector("list", nobs)
  M_suml <- vector("list", nobs)
  M_sumll <- vector("list", nobs)
  
  for (t in 1:nobs) {
    K_t <- kernelt(t, nobs, kernel_type = kernel_type, bandwidth = bw)
    
    M_cross <- matrix(0, nrow = ncol(Zmat), ncol = ncol(Zmat))
    
    for (j in 1:(2 * bw)) {
      w_j <- bartlett_kernel(j, L = lag)
      
      s_start <- max(1, t - bw)
      s_end   <- min(t + bw - j, nobs - j)
      
      if (s_end >= s_start) {
        for (s in s_start:s_end) {
          m_s  <- as.matrix(Zmat[s, ]) *
            (Y[s] - as.numeric(Xmat[s, ] %*% beta_mat[, s]))
          m_sj <- as.matrix(Zmat[s + j, ]) *
            (Y[s + j] - as.numeric(Xmat[s + j, ] %*% beta_mat[, s + j]))
          
          M_cross <- M_cross +
            w_j * diag(K_t)[s] * diag(K_t)[s + j] *
            (m_s %*% t(m_sj))
        }
      }
    }
    
    M_crossl[[t]] <- 2 * M_cross
    
    M_sumll[[t]] <- (M2l[[t]] + M_crossl[[t]])
    M_suml[[t]] <- solve(M_sumll[[t]])
    
    Q_l[[t]] <- (Ht[[t]] %*% M_sumll[[t]] %*% t(Ht[[t]]))
  }
  
  return(list(
    beta = beta_mat,
    var  = Q_l,
    opt_gmmw = M_suml,
    bw = bw
  ))
}


# ========== Example Usage ==========
set.seed(123)
T <- 100
data <- dgp_e2(T)  

X <- as.matrix(data$X_t)
Y <- as.matrix(data$Y_t)
W <- as.matrix(data[, c("Z_t1", "Z_t2", "Z_t3", "Z_t4")])
Z <- matrix(0, nrow = T, ncol = 0)


stage1 <- lgmm(
  X = X, Y = Y, Z = Z, W = W, bw = NULL,
  gmmw = NULL, withZ = 0,
  lag = 0, kernel_type = "epanechnikov"
)

gmw <- stage1$opt_gmmw
bw2 <- stage1$bw

stage2 <- lgmm(
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
