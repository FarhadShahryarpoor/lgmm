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





mirror_data <- function(data, bw) {
  if (is.vector(data)) {
    n <- length(data)
    top <- data[seq(bw, 1)]
    bottom <- data[seq(n, n - bw + 1)]
    return(c(top, data, bottom))
  } else if (is.matrix(data)) {
    n <- nrow(data)
    top <- data[seq(bw, 1), , drop = FALSE]
    bottom <- data[seq(n, n - bw + 1), , drop = FALSE]
    return(rbind(top, data, bottom))
  } else {
    stop("Data must be either a vector or a matrix.")
  }
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
  bw_list <- seq(0.3, 0.95, by = 0.02)
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
    print(t)
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
      
      if ((s_end >= s_start) & (lag != 0)) {
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



Optbw.bclgmm <- function(Xmat, Y, Zmat, gmmwl, kernel_type = "epanechnikov") {
  bw_list <- seq(0.4, 0.95, by = 0.02)
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
      if (s_end >= s_start) {
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








