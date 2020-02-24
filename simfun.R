########################################
## PURPOSE: Simulation Functions      ##
## BY: Kevin Josey                    ##
########################################

library(cbal)

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

  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  
  # effect coefficients
  beta <- c(210, 27.4, 13.7, 13.7, 13.7)
  gamma <- c(20, -13.7, 0, 0, 13.7)
  lambda <- c(0, -1, 0.5, -0.25, -0.1)
  theta <- c(0, 1, 0.3, -0.4, 0.1)
  
  # propensity score
  if (z_scen == "b") {
    e_X <- 1/(1 + exp( -( U %*% lambda) ) )
  } else { # z_scen == "a"
    e_X <- 1/(1 + exp( -( X %*% lambda) ) )
  }
  
  # Trial Participation
  if (s_scen == "b") {
    f_X <- 1/(1 + exp( -( U %*% theta) ) )
  } else { # s_scen == "a"
    f_X <- 1/(1 + exp( -( X %*% theta) ) )
  }  
  
  z <- rbinom(n, 1, e_X)
  s <- rbinom(n, 1, f_X)
  
  if (y_scen == "b") {
    W <- U
  } else { # y_scen == "b"
    W <- X
  }
  
  # outcome mean
  mu_0 <- W %*% beta
  mu_1 <- W %*% (beta + gamma)
  
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
  
  # Calibration - Transport
  
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
  
  target <- n_1*colMeans(X0)
  fit_tmp_t <- cfit(cmat = X1, target = target, distance = "entropy")
  
  q <- fit_tmp_t$weights
  cmat <- cbind(as.matrix((2*Z1 - 1)*X1), as.matrix( Z1*X1 ))
  b <- c(rep(0, times = ncol(X1)), c(t(X1) %*% q))
  fit_t <- try( cfit_transport(cmat = cmat, target = b, base_weights = q), silent = TRUE )
  
  if (!inherits(fit_t, "try-error")) {
    
    A <- cbind((2*Z - 1)*X, Z*X)
    p <- rep(1, n_1 + n_0)
    p[S == 1] <- fit_t$weights
    q <- rep(1, n_1 + n_0)
    q[S == 1] <- fit_t$base_weights
    coefs_1 <- fit_t$coefs
    coefs_2 <- fit_tmp_t$coefs
    
    tau_t <- sum(S*(2*Z - 1)*p*q*Y)/sum(S*p*q*Z)
    
    dp <- as.vector( -exp(-A %*% coefs_1) )
    dq <- as.vector( -exp(-X %*% coefs_2) )
    U <- matrix(0, ncol = 4*m, nrow = 4*m)
    v <- rep(0, times = 4*m + 1)
    meat <- matrix(0, ncol = 4*m + 1, nrow = 4*m + 1)
    
    for (i in 1:(n_1 + n_0)) {
      
      U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] + 
        S[i]* q[i]*dp[i] * (A[i,] %*% t(A[i,]))
      U[1:m,(2*m + 1):(3*m)] <- U[1:m,(2*m + 1):(3*m)] + 
        S[i]*(2*Z[i] - 1)* dq[i]*p[i] * (X[i,] %*% t(X[i,]))
      U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] +
        S[i]*Z[i]* dq[i]*p[i] * (X[i,] %*% t(X[i,])) - dq[i] * (X[i,] %*% t(X[i,]))
      U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] + 
        S[i]*dq[i]*X[i,] %*% t(X[i,])
      
      U[(2*m + 1):(3*m),(3*m + 1):(4*m)] <- U[(2*m + 1):(3*m),(3*m + 1):(4*m)] + S[i]*diag(-1, m, m)
      U[(3*m + 1):(4*m),(3*m + 1):(4*m)] <- U[(3*m + 1):(4*m),(3*m + 1):(4*m)] + (1 - S[i])*diag(-1, m, m)
      
      v[1:(2*m)] <- v[1:(2*m)] + S[i]*(2*Z[i] - 1)*q[i]*dp[i] * (Y[i] - Z[i]*tau_t) * A[i,]
      v[(2*m+1):(3*m)] <- v[(2*m + 1):(3*m)] + S[i]*(2*Z[i] - 1) * dq[i]*p[i] * (Y[i] - Z[i]*tau_t) * X[i,]
      v[4*m + 1] <- v[4*m + 1] - S[i]*q[i]*p[i]*Z[i]
      meat <- meat + tcrossprod(esteq_transport(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], weights = p[i],
                                                base_weights = q[i], target = target/n_1, tau = tau_t))
      
      
    }
    
    invbread <- matrix(0, nrow = 4*m + 1, ncol = 4*m + 1)
    invbread[1:(4*m),1:(4*m)] <- U
    invbread[4*m + 1,] <- v
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error"))
      bread <- try(MASS::ginv(invbread))
    
    if (inherits(bread, "try-error")) {
      
      var_t <- NA
      cp_t <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread) 
      var_t <- sandwich[4*m + 1, 4*m + 1]
      cp_t <- as.numeric(tau_t - sqrt(var_t)*1.96 <= PATE & tau_t + sqrt(var_t)*1.96 >= PATE)
      
    }
    
  } else {
    
    tau_t <- NA
    var_t <- NA
    cp_t <- NA
    
  }
  
  # Calibration - Fusion
  
  target <- n_1*colMeans(X0)
  fit_tmp_f <- cfit(cmat = S*X, target = target, distance = "entropy")

  q <- fit_tmp_f$weights
  # cmat <- cbind(as.matrix((2*Z - 1)*X), as.matrix( Z*X ))
  # b <- c(rep(0, times = ncol(X)), c(t(X) %*% q))
  cmat <- cbind(as.matrix((2*Z - 1)*S*X), as.matrix( Z*S*X ), as.matrix((2*Z - 1)*(1 - S)*X), as.matrix( Z*(1 - S)*X ))
  b <- c(rep(0, times = ncol(X1)), c(t(X1) %*% q[S == 1]), rep(0, times = ncol(X0)), c(t(X0) %*% q[S == 0]))
  fit_f <- try( cfit_transport(cmat = cmat, target = b, base_weights = q), silent = TRUE )

  if (!inherits(fit_f, "try-error")) {

    A <- cbind((2*Z - 1)*X, Z*X)
    p <- fit_f$weights
    q <- fit_f$base_weights
    coefs_1 <- fit_f$coefs
    coefs_2 <- fit_tmp_f$coefs

    cdat <- as.data.frame(cbind(Y0 = Y0, unname(X0[,-1])))
    m_0 <- predict(lm(Y0 ~ ., data = cdat[Z0 == 0,]), newdata = as.data.frame(unname(X)))
    m_1 <- predict(lm(Y0 ~ ., data = cdat[Z0 == 1,]), newdata = as.data.frame(unname(X)))

    tau_f <- sum(p*q*(2*Z - 1)*Y)/sum(p*q*Z)
    var_f <- NA
    cp_f <- NA

  #   dp <- as.vector( -exp(-A %*% coefs_1) )
  #   dq <- as.vector( -exp(-X %*% coefs_2) )
  #   U <- matrix(0, ncol = 4*m, nrow = 4*m)
  #   v <- rep(0, times = 4*m + 1)
  #   meat <- matrix(0, ncol = 4*m + 1, nrow = 4*m + 1)
  # 
  #   for (i in 1:(n_1 + n_0)) {
  # 
  #     U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] + S[i]*q[i]*dp[i] * (A[i,] %*% t(A[i,]))
  #     U[1:m,(2*m + 1):(3*m)] <- U[1:m,(2*m + 1):(3*m)] + S[i]*(2*Z[i] - 1)*dq[i]*p[i] * (X[i,] %*% t(X[i,]))
  #     U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] + S[i]*Z[i]*dq[i]*p[i] * (X[i,] %*% t(X[i,])) - dq[i] * (X[i,] %*% t(X[i,]))
  #     U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] + S[i]*dq[i]*X[i,] %*% t(X[i,])
  # 
  #     U[(2*m + 1):(3*m),(3*m + 1):(4*m)] <- U[(2*m + 1):(3*m),(3*m + 1):(4*m)] + S[i]*diag(-1, m, m)
  #     U[(3*m + 1):(4*m),(3*m + 1):(4*m)] <- U[(3*m + 1):(4*m),(3*m + 1):(4*m)] + (1 - S[i])*diag(-1, m, m)
  # 
  #     v[1:(2*m)] <- v[1:(2*m)] + S[i]*(2*Z[i] - 1)*q[i]*dp[i] * (Y[i] - Z[i]*tau_f) * A[i,]
  #     v[(2*m+1):(3*m)] <- v[(2*m + 1):(3*m)] + S[i]*(2*Z[i] - 1) * dq[i]*p[i] * (Y[i] - Z[i]*tau_f) * X[i,]
  #     v[4*m + 1] <- v[4*m + 1] - S[i]*q[i]*p[i]*Z[i]
  #     meat <- meat + tcrossprod(esteq_transport(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], weights = p[i],
  #                                               base_weights = q[i], target = target/n_1, tau = tau_f))
  # 
  # 
  #   }
  # 
  #   invbread <- matrix(0, nrow = 4*m + 1, ncol = 4*m + 1)
  #   invbread[1:(4*m),1:(4*m)] <- U
  #   invbread[4*m + 1,] <- v
  # 
  #   bread <- try(solve(invbread), silent = TRUE)
  # 
  #   if (inherits(bread, "try-error"))
  #     bread <- try(MASS::ginv(invbread))
  # 
  #   if (inherits(bread, "try-error")) {
  # 
  #     var_f <- NA
  #     cp_f <- NA
  # 
  #   } else {
  # 
  #     sandwich <- bread %*% meat %*% t(bread)
  #     var_f <- sandwich[4*m + 1, 4*m + 1]
  #     cp_f <- as.numeric(tau_f - sqrt(var_f)*1.96 <= PATE & tau_f + sqrt(var_f)*1.96 >= PATE)
  # 
  #   }
    
  } else {

    tau_f <- NA
    var_f <- NA
    cp_f <- NA

  }

  # IOSW
  
  glmdat <- data.frame(S = S, X = X[,-1])
  glmfit <- glm(S ~ ., data = glmdat, family = binomial(link = "logit"))
  glmprob <- glmfit$fitted.values
  glmwts_tmp <- (1 - glmprob)/(glmprob)
  glmwts <- glmwts_tmp[S == 1]
  
  baldat <- data.frame(Z = Z1, X = X1[,-1])
  balfit <- glm(Z ~ ., data = baldat, family = binomial(link = "logit"))
  balprob <- balfit$fitted.values
  balwts <- ifelse(Z1 == 1, 1/balprob, 1/(1 - balprob))
  
  wts <- ifelse(Z1 == 1, glmwts/balprob, glmwts/(1 - balprob))
  
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
  
  # GLM Results
  
  design <- svydesign(ids = ~ 1, weights = ~ wts, data = data.frame(wts = wts, Y1 = Y1, Z1 = Z1))
  smod <- svyglm(Y1 ~ Z1, design = design, family = gaussian)
  glmrslt <- coef(smod)[2]
  
  # Outcome Results
  
  outrslt <- mean(g_1 - g_0)
  
  # TMLE Results
  
  tmlerslt <- tmleest$estimate
  
  # Combine Results
  
  tau <- c(glmrslt, outrslt, tmlerslt, tau_t. tau_f)
  cp <- c(cp_t)
  
  return(list(tau = tau, cp = cp, PATE = PATE))
  
}
