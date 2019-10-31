########################################
## PURPOSE: Simulation Functions Code ##
## BY: Kevin Josey                    ##
########################################

gen_data <- function(n_0, n_1, prob, sig2, rho = 0, scenario = c("baseline", "positivity", "interaction", "sparse")) {
  
  # error variance
  R <- matrix(rho, nrow = 2, ncol = 2)
  diag(R) <- 1
  V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
  Sig <- V %*% R %*% V
  
  # treatment assignment
  z0 <- rbinom(n_0, 1, prob)
  z1 <- rbinom(n_1, 1, prob)
  
  if (scenario == "baseline") {
    
    # effect coefficients
    beta <- c(10, -1, -2, 2, 1)
    # beta <- c(10, 1, 2, 1)
    alpha <- c(2, 2, 1, -1, -2)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, size = 1, prob = 0.5)
    
    # covariate values
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, size = 1, prob = 0.5)
   
    theta <- c(1, -1, 0.6, 0, 0.5)
    X0 <- cbind(int0, x00, x01, x02, x03)
    X1 <- cbind(int1, x10, x11, x12, x13)
    
    # Outcome Mean
    mu_10 <- X1 %*% beta
    mu_11 <- X1 %*% beta + X1 %*% alpha
    
    SATE <- mean(X0 %*% alpha)
    PATE <- as.numeric(t(theta) %*% alpha)
    
  } else if (scenario == "positivity") {
    
    # effect coefficients
    beta <- c(10, -1, -2, 2, 1)
    alpha <- c(2, 2, 1, -1, -2)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 1)
    x01 <- rbinom(n_0, size = 1, prob = 0.7)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, size = 1, prob = 0.5)
    
    # covariate values
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.3)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, size = 1, prob = 0.5)
    
    theta <- c(1, -1, 0.7, 0, 0.5)
    X0 <- cbind(int0, x00, x01, x02, x03)
    X1 <- cbind(int1, x10, x11, x12, x13)
    
    # Outcome Mean
    mu_10 <- X1 %*% beta
    mu_11 <- X1 %*% beta + z1*X1 %*% alpha
    
    SATE <- mean(X0 %*% alpha)
    PATE <- as.numeric(t(theta) %*% alpha)
    
  } else if (scenario == "interaction") {
    
    # effect coefficients
    beta <- c(10, -1, -2, 2, 1, 2, -1)
    alpha <- c(2, 2, 1, -1, -2, -1, 2)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, size = 1, prob = 0.5)
    
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, size = 1, prob = 0.5)
    
    theta <- c(1, -1, 0.6, 0, 0.5, -0.6, -0.5)
    X0_tmp <- cbind(int0, x00, x01, x02, x03, x00*x01, x00*x03)
    X0 <- cbind(int0, x00, x01, x02, x03)
    X1_tmp <- cbind(int1, x10, x11, x12, x13, x10*x11, x10*x13)
    X1 <- cbind(int1, x10, x11, x12, x13)
    
    mu_10 <- X1_tmp %*% beta
    mu_11 <- X1_tmp %*% beta + X1_tmp %*% alpha
    
    SATE <- mean(X0_tmp %*% alpha)
    PATE <- as.numeric(t(theta) %*% alpha)
    
  } else if (scenario == "sparse") {
    
    # effect coefficients
    
    beta <- c(10, -1, -2, 2, 1, rep(0, times = 4))
    alpha <- c(2, 2, 1, -1, -2, rep(0, times = 4))
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, 1, prob = 0.5)

    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, 1, prob = 0.5)
    
    x04 <- rnorm(n_0, mean = 1, sd = 2)
    x05 <- rbinom(n_0, size = 1, prob = 0.4)
    x06 <- rnorm(n_0, mean = 0, sd = 1)
    x07 <- rbinom(n_0, size = 1, prob = 0.5)

    x14 <- rnorm(n_1, mean = -1, sd = 2)
    x15 <- rbinom(n_1, size = 1, prob = 0.6)
    x16 <- rnorm(n_1, mean = 0, sd = 1)
    x17 <- rbinom(n_1, size = 1, prob = 0.5)
    
    theta <- c(1, -1, 0.6, 0, 0.5, 1, 0.4, 0, 0.5)
    X0 <- cbind(int0, x00, x01, x02, x03, x04, x05, x06, x07)
    X1 <- cbind(int1, x10, x11, x12, x13, x14, x15, x16, x17)
    
    mu_10 <- X1%*%beta
    mu_11 <- X1%*%beta + X1%*%alpha
    
    SATE <- mean(X0 %*% alpha)
    PATE <- as.numeric(t(theta) %*% alpha)
    
  }

  # potential outcomes
  eval <- eigen(Sig, symmetric = TRUE)
  y1_init <- matrix(stats::rnorm(n_1*2, 0, 1), nrow = n_1, ncol = 2) # iid potential outcomes
  y1_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y1_init)) # SVD
  y1_pot <- y1_tmp + cbind(mu_10, mu_11) # include causal effect
  
  # observed outcome
  y1 <- z1*y1_pot[,2] + (1 - z1)*y1_pot[,1]

  out <- list(X0 = X0, X1 = X1, y1 = y1, z0 = z0, z1 = z1, SATE = SATE, PATE = PATE)
  return(out)
  
}

# Fits the balancing weights using a variety of methods
simfit <- function(idx = 1, simDat) {
  
  dat <- simDat[,idx]
  SATE <- dat$SATE
  PATE <- dat$PATE
  n_1 <- nrow(dat$X1)
  
  # Entropy Calibration
  
  tm <- colMeans(dat$X0)
  cm <- dat$X1
  
  entfit_1 <- calib(X = cm, Z = dat$z1, target = tm, complex = TRUE)
  entest_1 <- try( calest(obj = entfit_1, Y = dat$y1, method = "sandwich"), silent = TRUE )
  
  entfit_2 <- calib(X = cm, Z = dat$z1, target = tm, complex = FALSE)
  entest_2 <- try( calest(obj = entfit_2, Y = dat$y1, method = "sandwich"), silent = TRUE )
  
  # IOSW
  
  baldat <- rbind(cbind(z = dat$z0, dat$X0[,-1]), cbind(z = dat$z1, dat$X1[,-1]))
  glmdat <- data.frame(S = c(rep(0, times = nrow(dat$X0)), rep(1, times = nrow(dat$X1))), baldat)
  glmfit <- glm(S ~ ., data = glmdat, family = binomial(link = "logit"))
  
  glmprob <- glmfit$fitted.values
  glmwts <- (1 - glmprob)/(glmprob)
  wts <- glmwts[glmdat$S == 1]/(dat$z1*mean(dat$z1) + (1 - dat$z1)*mean(1 - dat$z1))
  
  # OM
  cdat <- as.data.frame(cbind(y1 = dat$y1, unname(dat$X1[,-1])))
  g_0 <- predict(lm(y1 ~ ., data = cdat[dat$z1 == 0,]), newdata = as.data.frame(unname(dat$X0)))
  g_1 <- predict(lm(y1 ~ ., data = cdat[dat$z1 == 1,]), newdata = as.data.frame(unname(dat$X0)))
  
  # AIPW
  mod_1 <- lm(y1 ~ ., data = cdat[dat$z1 == 1,])
  m_1 <- predict(mod_1, cdat)
  mod_0 <- lm(y1 ~ ., data = cdat[dat$z1 == 0,])
  m_0 <- predict(mod_0, cdat)        
  AIPW_1 <- sum(dat$z1*wts*(dat$y1 - m_1))
  AIPW_0 <- sum((1-dat$z1)*wts*(dat$y1 - m_0))
  norm_1 <- (sum(dat$z1*wts))^-1
  norm_0 <- (sum((1-dat$z1)*wts))^-1

  # Entropy Results
  
  if (inherits(entest_1, "try-error")){
    
    entrslt <- NA
    entvar <- NA
    entcps <- NA
    entcpp <- NA
    
  } else { 
    
    entrslt <- entest_1$tau
    entvar <- entest_1$variance
    entcps <- as.numeric(entrslt - sqrt(entvar)*1.96 <= SATE & entrslt + sqrt(entvar)*1.96 >= SATE)
    entcpp <- as.numeric(entrslt - sqrt(entvar)*1.96 <= PATE & entrslt + sqrt(entvar)*1.96 >= PATE)
    
  }
  
  # MOM Result
  
  if (inherits(entest_2, "try-error"))
    momrslt <- NA
  else
    momrslt <- entest_2$tau
  
  
  # GLM Results
  
  design <- svydesign(ids = ~ 1, weights = ~ wts, data = data.frame(wts = wts, Y1 = dat$y1, Z1 = dat$z1))
  smod <- svyglm(Y1 ~ Z1, design = design, family = gaussian)
  glmrslt <- coef(smod)[2]
  
  # Outcome Results
  
  outrslt <- mean(g_1 - g_0)
  outvar <- as.numeric(t(tm) %*% vcov(lm(y1 ~ ., data = cdat[dat$z1 == 0,])) %*% tm + 
    t(tm) %*% vcov(lm(y1 ~ ., data = cdat[dat$z1 == 1,])) %*% tm)
  outcps <- as.numeric(outrslt - sqrt(outvar)*1.96 <= SATE & outrslt + sqrt(outvar)*1.96 >= SATE)
  outcpp <- as.numeric(outrslt - sqrt(outvar)*1.96 <= PATE & outrslt + sqrt(outvar)*1.96 >= PATE)
  
  # AIPW Results
  
  arg_1 <- norm_1*AIPW_1 + mean(g_1)
  arg_0 <- norm_0*AIPW_0 + mean(g_0)
  aipwrslt <- arg_1 - arg_0
  
  # Combine Results
  
  tau <- c(glmrslt, outrslt, aipwrslt, momrslt, entrslt)
  cp <- c(entcps, entcpp, outcps, outcpp)
  
  return(list(tau = tau, cp = cp))
  
}
