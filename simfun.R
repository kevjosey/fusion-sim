########################################
## PURPOSE: Simulation Functions      ##
## BY: Kevin Josey                    ##
########################################

hte_data <- function(n, sig2 = 5, y_scen = c("a", "b"), z_scen = c("a", "b"), s_scen = c("a", "b")){
  
  # error variance
  R <- matrix(0, nrow = 2, ncol = 2)
  diag(R) <- 1
  V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
  Sig <- V %*% R %*% V
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  
  # transformed predictors
  u1 <- as.numeric(scale(exp(x1/2)))
  u2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
  u3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
  u4 <- as.numeric(scale((x2 + x4 + 20)^2))
  
  v1 <- stats::rnorm(n, 0, 1)
  v2 <- stats::rnorm(n, 0, 1)
  v3 <- stats::rnorm(n, 0, 1)
  v4 <- stats::rnorm(n, 0, 1)
  
  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  V <- cbind(int = rep(1, n), v1, v2, v3, v4)
  
  # effect coefficients
  beta <- c(210, 27.4, 13.7, 13.7, 13.7)
  alpha <- c(20, -13.7, 0, 0, 13.7)
  lambda <- c(0.1, -1, 0.5, -0.25, -0.1)
  delta <- c(0, 1, 0.5, 0.4, -0.2)
  gamma <- c(0, 1, 0.3, -0.4, 0.2)
  
  # Trial Participation
  if (s_scen == "b") {
    f_X <- 1/(1 + exp( -( U %*% gamma) ) )
  } else { # s_scen == "a"
    f_X <- 1/(1 + exp( -( X %*% gamma) ) )
  } 
  
  s <- rbinom(n, 1, f_X)
  
  # propensity score
  if (z_scen == "b") {
    e_X <- 1/(1 + exp( -( U %*% lambda) ) )
  } else { # z_scen == "a"
    e_X <- s/(1 + exp( -( U %*% lambda) ) ) + (1 - s)/(1 + exp( -( X %*% delta) ) )
  }
  
  z <- rbinom(n, 1, e_X)
  
  if (y_scen == "b") {
    W <- U
  } else { # y_scen == "b"
    W <- X
  }
  
  # outcome mean
  mu_0 <- W %*% beta
  mu_1 <- W %*% (alpha + beta)
  
  # potential outcomes
  eval <- eigen(Sig, symmetric = TRUE)
  y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
  y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
  y_pot <- y_tmp + cbind(mu_0, mu_1) # include causal effect
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  PATE <- mean(y_pot[s == 0,2] - y_pot[s == 0,1])
  
  # create simulation dataset
  sim <- list(y = y, z = z, X = X, U = U, s = s, PATE = PATE)
  
  return(sim)
  
}

# Fits the balancing weights using a variety of methods
simfit <- function(idx = 1, simDat) {
  
  dat <- simDat[,idx]
  PATE <- dat$PATE
  
  # stuff
  
  S <- dat$s
  Y <- dat$y
  Y1 <- dat$y[S == 1]
  Y0 <- dat$y[S == 0]
  Z <- dat$z
  Z1 <- dat$z[S == 1]
  Z0 <- dat$z[S == 0]
  X <- dat$X
  X1 <- dat$X[S == 1,]
  X0 <- dat$X[S == 0,]
  
  n_1 <- sum(S)
  n_0 <- sum(1 - S)
  m <- ncol(X)
  
  tm <- colMeans(X0)
  
  R <- as.matrix(S*X)
  fit_tmp <- cfit(cmat = R, target = n_1*tm, distance = "entropy")
  q <- fit_tmp$weights
  
  # Calibration - Transport
  
  A <- cbind(as.matrix(S*(2*Z - 1)*X), as.matrix( S*Z*X ))
  b <- c(rep(0, times = ncol(X)), c(t(R) %*% q))
  fit_t <- try( cfit(cmat = A, target = b, base_weights = q, distance = "transport"), silent = TRUE )
  
  if (!inherits(fit_t, "try-error")) {
    
    p <- fit_t$weights
    coefs_1 <- fit_t$coefs
    coefs_2 <- fit_tmp$coefs
    
    tau_t <- sum(S*(2*Z - 1)*p*q*Y)/sum(S*p*q*Z)
    
    dp <- as.vector( -exp(-A %*% coefs_1) )
    dq <- as.vector( -exp(-R %*% coefs_2) )
    U <- matrix(0, ncol = 3*m, nrow = 3*m)
    v <- rep(0, times = 3*m + 1)
    meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)
    
    for (i in 1:(n_1 + n_0)) {
      
      U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] + 
        q[i]*dp[i] * (A[i,] %*% t(A[i,]))
      
      U[1:m,(2*m + 1):(3*m)] <- U[1:m,(2*m + 1):(3*m)] + 
        dq[i]*p[i] * (2*Z[i] - 1) * (X[i,] %*% t(R[i,]))
      U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] +
        dq[i]*(p[i]*Z[i]*X[i,] - X[i,]) %*% t(R[i,])
      U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] + 
        dq[i] * (X[i,] %*% t(R[i,]))
      
      v[1:(2*m)] <- v[1:(2*m)] + q[i]*dp[i] * (2*Z[i] - 1)*(Y[i] - Z[i]*tau_t)*A[i,]
      v[(2*m+1):(3*m)] <- v[(2*m + 1):(3*m)] + dq[i]*p[i] * (2*Z[i] - 1)*(Y[i] - Z[i]*tau_t)*R[i,]
      v[3*m + 1] <- v[3*m + 1] - S[i]*q[i]*p[i]*Z[i]
      meat <- meat + tcrossprod(esteq_transport(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                                weights = p[i], base_weights = q[i], tau = tau_t,
                                                n_1 = n_1, n_0 = n_0))
      
      
    }
    
    invbread <- matrix(0, nrow = 3*m + 1, ncol = 3*m + 1)
    invbread[1:(3*m),1:(3*m)] <- U
    invbread[3*m + 1,] <- v
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error"))
      var_t <- NA
    
    else {
      
      sandwich <- bread %*% meat %*% t(bread)
      var_t <- sandwich[3*m + 1, 3*m + 1]
      
    }
    
  } else {
    
    tau_t <- NA
    var_t <- NA
    
  }
  
  # Calibration - Fusion

  A <- cbind(as.matrix((2*Z - 1)*X), as.matrix( Z*X ))
  b <- c(rep(0, times = ncol(X)), c(t(X) %*% q))
  fit_f <- try( cfit(cmat = A, target = b, base_weights = q, distance = "transport"), silent = TRUE )

  if (!inherits(fit_f, "try-error")) {
    
    p <- fit_f$weights
    coefs_1 <- fit_f$coefs
    coefs_2 <- fit_tmp$coefs

    tau_f <- sum(p*q*(2*Z - 1)*Y)/sum(p*q*Z)

    dp <- as.vector( -exp(-A %*% coefs_1) )
    dq <- as.vector( -exp(-R %*% coefs_2) )
    U <- matrix(0, ncol = 3*m, nrow = 3*m)
    v <- rep(0, times = 3*m + 1)
    meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)

    for (i in 1:(n_1 + n_0)) {

      U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] + q[i]*dp[i] * (A[i,] %*% t(A[i,]))
      
      U[1:m,(2*m + 1):(3*m)] <- U[1:m,(2*m + 1):(3*m)] + 
        dq[i]*p[i] * (2*Z[i] - 1) * (X[i,] %*% t(R[i,]))
      U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] + 
        dq[i]*(p[i]*Z[i]*X[i,] - X[i,]) %*% t(R[i,])
      U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] + 
        dq[i]*(X[i,] %*% t(R[i,]))

      v[1:(2*m)] <- v[1:(2*m)] + q[i]*dp[i] * (2*Z[i] - 1)*(Y[i] - Z[i]*tau_f)*A[i,]
      v[(2*m+1):(3*m)] <- v[(2*m + 1):(3*m)] + dq[i]*p[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau_f)*R[i,]
      v[3*m + 1] <- v[3*m + 1] - q[i]*p[i]*Z[i]
      meat <- meat + tcrossprod(esteq_fusion(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                             weights = p[i], base_weights = q[i], tau = tau_f,
                                             n_1 = n_1, n_0 = n_0))


    }

    invbread <- matrix(0, nrow = 3*m + 1, ncol = 3*m + 1)
    invbread[1:(3*m),1:(3*m)] <- U
    invbread[3*m + 1,] <- v

    bread <- try(solve(invbread), silent = TRUE)

    if (inherits(bread, "try-error"))
      var_f <- NA
    
    else {

      sandwich <- bread %*% meat %*% t(bread)
      var_f <- sandwich[3*m + 1, 3*m + 1]

    }
    
  } else {

    tau_f <- NA
    var_f <- NA

  }
  
  # OM
  
  cdat <- as.data.frame(cbind(Y1 = Y1, unname(X1[,-1])))
  g_0 <- predict(lm(Y1 ~ ., data = cdat[Z1 == 0,]), newdata = as.data.frame(unname(X0)))
  g_1 <- predict(lm(Y1 ~ ., data = cdat[Z1 == 1,]), newdata = as.data.frame(unname(X0)))
  
  # TMLE
  
  smod <- "S ~ x1 + x2 + x3 + x4" 
  ymod <- "Y ~ x1 + x2 + x3 + x4 + Z*x1 + Z*x2 + Z*x3 + Z*x4" 
  zmod <- "Z ~ x1 + x2 + x3 + x4" 
  
  W <- data.frame(X[,-1], S = S, Y = Y, Z = Z)
  
  tmleest <- try( tmle(S = S, Y = Y, Z = Z, data = as.data.frame(W), nsitemodel = smod, nzmodel = zmod, noutmodel = ymod), silent = TRUE )
  
  # Augmented Calibration
  
  cdat_2 <- as.data.frame(cbind(Y0 = Y0, unname(X0[,-1])))
  m_0 <- predict(lm(Y0 ~ ., data = cdat_2[Z0 == 0,]), newdata = as.data.frame(unname(X)))
  m_1 <- predict(lm(Y0 ~ ., data = cdat_2[Z0 == 1,]), newdata = as.data.frame(unname(X)))
  
  glmdat <- data.frame(S = S, X = X[,-1])
  glmfit <- glm(S ~ ., data = glmdat, family = binomial(link = "logit"))
  glmprob <- glmfit$fitted.values
  glmwts_tmp <- (1 - glmprob)/(glmprob)
  glmwts <- glmwts_tmp[S == 1]
  
  baldat <- data.frame(Z = Z1, X = X1[,-1])
  balfit <- glm(Z ~ ., data = baldat, family = binomial(link = "logit"))
  balprob <- balfit$fitted.values
  balwts <- ifelse(Z1 == 1, 1/balprob, 1/(1 - balprob))
  ipw <- rep(1, times = n)
  ipw[S == 1] <- balwts
  
  tau_e <- sum(S*q*(Z*(Y - m_1)/ipw - (1 - Z)*(Y - m_0)/ipw))/n_1 + sum((1-S)*(m_1 - m_0))/n_0
  
  # Outcome Results
  
  outrslt <- mean(g_1 - g_0)
  
  # TMLE Results
  
  tmlerslt <- tmleest$estimate
  tmlevar <- tmleest$variance
  
  # Combine Results
  
  tau <- c(outrslt, tmlerslt, tau_t, tau_e, tau_f)
  se <- c(sqrt(var_t), sqrt(var_f))
  
  return(list(tau = tau, se = se, PATE = PATE))
  
}
