########################################
## PURPOSE: Simulation Functions      ##
## BY: Kevin Josey                    ##
########################################

gen_data <- function(n, sig2 = 5, y_scen = c("a", "b"), z_scen = c("a", "b"), s_scen = c("a", "b")){
  
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
  u1 <- as.numeric(scale(exp((x1 + x4))))
  u2 <- as.numeric(scale((x1 + x2)^3))
  u3 <- as.numeric(scale((x2 + x3)^2))
  u4 <- as.numeric(scale(log(abs(x3*x4))))
  
  # u1 <- as.numeric(scale(exp(x1/2)))
  # u2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
  # u3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
  # u4 <- as.numeric(scale((x2 + x4 + 20)^2))
  
  # create matrix
  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  # V <- cbind(int = rep(1, n), v1, v2, v3, v4)
  
  # coefficients
  beta <- c(5, -1, 3, -3, 1)
  alpha <- c(10, -3, -1, 1, 3)
  lambda <- c(0, 0, 1, -0.5, 0.5)
  delta <- c(0, -1, -0.5, 0, 0.5)
  gamma <- c(0, -0.5, -1, -0.5, 1)
  
  # Trial Participation
  if (s_scen == "b") {
    f_X <- 1/(1 + exp( -( U %*% gamma) ) )
  } else { # s_scen == "a"
    f_X <- 1/(1 + exp( -( X %*% gamma) ) )
  } 
  
  s <- rbinom(n, 1, f_X)
  
  # propensity score
  if (z_scen == "b") {
    e_X <- s/(1 + exp( -( U %*% delta) ) ) + (1 - s)/(1 + exp( -( U %*% lambda) ) )
  } else { # z_scen == "a"
    e_X <- s/(1 + exp( -( X %*% delta) ) ) + (1 - s)/(1 + exp( -( X %*% lambda) ) )
  }
  
  z <- rbinom(n, 1, e_X)
  
  if (y_scen == "b") {
    W <- U
  } else { # y_scen == "b"
    W <- X
  }
  
  # outcome mean
  mu_0 <- W %*% beta
  mu_1 <- W %*% alpha
  
  # potential outcomes
  eval <- eigen(Sig, symmetric = TRUE)
  y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
  y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
  y_pot <- y_tmp + cbind(mu_0, mu_1) # include causal effect
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  PATE <- mean(y_pot[s == 0,2] - y_pot[s == 0,1])
  
  # create simulation dataset
  sim <- list(y = y, z = z, X = X, s = s, PATE = PATE)
  
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
  n <- n_1 + n_0
  m <- ncol(X)
  
  theta <- colMeans(X0)
  
  # Calibration - Transport
  
  A <- cbind(as.matrix(S*Z*X), as.matrix( S*(1 - Z)*X ))
  b <- c(n_1*theta, n_1*theta)
  fit_t <- cfit(cmat = A, target = b, distance = "entropy")
  est_t <- transport_estimate(obj = fit_t, S = S, X = X, Z1 = Z1, Y1 = Y1)
  cal_t <- est_t$estimate
  var_t <- est_t$variance
  
  # Calibration - Fusion
  
  R <- as.matrix(S*X)
  fit_base <- cfit(cmat = R, target = n_1*theta, distance = "entropy")
  q <- fit_base$weights
  
  A <- cbind(as.matrix(Z*X), as.matrix( (1 - Z)*X ), as.matrix(S*X))
  b <- c(n*theta, n*theta, n_1*theta)
  fit_f <- cfit(cmat = A, target = b, distance = "entropy")
  est_f <- fusion_estimate(obj = fit_f, base_obj = fit_base, S = S, X = X, Z = Z, Y = Y)
  cal_f <- est_f$estimate
  var_f <- est_f$variance
  
  # TMLE
  
  smod <- "S ~ x1 + x2 + x3 + x4" 
  ymod <- "Y ~ x1 + x2 + x3 + x4 + Z*x1 + Z*x2 + Z*x3 + Z*x4" 
  zmod <- "Z ~ x1 + x2 + x3 + x4" 
  
  W <- data.frame(X[,-1], S = S, Y = Y, Z = Z)
  
  tmleest_t <- try( tmle(S = S, Y = Y, Z = Z, data = as.data.frame(W), 
                       nsitemodel = smod, nzmodel = zmod, noutmodel = ymod, fusion = FALSE), silent = TRUE )
  
  tmleest_f <- try( tmle(S = S, Y = Y, Z = Z, data = as.data.frame(W), 
                       nsitemodel = smod, nzmodel = zmod, noutmodel = ymod, fusion = TRUE), silent = TRUE )
  
  # Augmented Calibration - Transport
  
  cdat <- as.data.frame(cbind(Y = Y, unname(X[,-1])))
  
  m_0 <- predict(lm(Y ~ ., data = cdat, subset = (S == 1 & Z == 0)), newdata = cdat)
  m_1 <- predict(lm(Y ~ ., data = cdat, subset = (S == 1 & Z == 1)), newdata = cdat)
  
  bdat <- data.frame(Z = Z, X = X[,-1])
  bprob <- predict(glm(Z ~ ., data = bdat, family = quasibinomial(link = "logit"), subset = S == 1), newdata = bdat, type = "response") 
  ipw <- ifelse(Z == 1, 1/bprob, 1/(1 - bprob))
  
  aug_t <- sum(S*q*(Z*(Y - m_1)/ipw - (1 - Z)*(Y - m_0)/ipw))/n_1 + sum((1-S)*(m_1 - m_0))/n_0
  
  # Augmented Calibration - Fusion
  
  g_0 <- predict(lm(Y ~ ., data = cdat, subset = Z == 0), newdata = cdat)
  g_1 <- predict(lm(Y ~ ., data = cdat, subset = Z == 1), newdata = cdat)
  
  aug_f <- sum(S*q*(Z*(Y - g_1)/ipw  - (1 - Z)*(Y - g_0)/ipw))/n_1 + sum((1 - S)*(g_1 - g_0))/n_0
  
  # TMLE Results
  
  tmle_t <- tmleest_t$estimate
  tmle_f <- tmleest_f$estimate
  
  # Combine Results
  
  tau <- c(tmle_t, aug_t, cal_t, aug_f, cal_f)
  se <- c(sqrt(var_t), sqrt(var_f))
  
  return(list(tau = tau, se = se, PATE = PATE))
  
}
