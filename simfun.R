########################################
## PURPOSE: Simulation Functions      ##
## BY: Kevin Josey                    ##
########################################

gen_data <- function(n, sig2 = 5, scenario = c("base", "exchange",
                                               "ps-mis", "out-mis", 
                                               "prognostic","instrument",
                                               "sample-overlap", "treat-overlap")){
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  
  # transformed predictors
  u1 <- as.numeric(scale(exp((x1 + x4)/2)))
  u2 <- as.numeric(scale(x2/(1 + exp(x1))))
  u3 <- as.numeric(scale(log(abs(x2*x3))))
  u4 <- as.numeric(scale((x3 + x4)^2))
  
  # v1 <- as.numeric(scale(exp(x1/2)))
  # v2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
  # v3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
  # v4 <- as.numeric(scale((x2 + x4 + 10)^2))
  
  # create matrix
  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  # V <- cbind(int = rep(1, n), v1, v2, v3, v4)
  
  # coefficients
  beta0 <- c(2, -3, -1, 1, 3)
  beta1 <- c(0, 2, -2, -2, 2)
  alpha <- c(-2, -1, 3, -3, 1)
  lambda <- c(0, 0.5, -0.5, 0.5, -0.5)
  delta <- c(-0.5, 0, 0, 0, 0)
  gamma <- c(0.5, -0.5, 0.5, -0.5, 0.5)
  
  # beta0 <- c(4, -3, -1, 1, 3)
  # beta1 <- c(2, -1, -3, 3, 1)
  # alpha <- c(2, 2, 2, -2, -2)
  # lambda <- c(0.5, -0.5, 0.5, -0.5, 0.5)
  # delta <- c(-0.25, 0, 0, 0, 0)
  # gamma <- c(-0.25, 0.25, -0.25, 0.25, -0.25)
  
  if (scenario == "base"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(1/(1 + exp(-X %*% lambda)))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(X %*% beta0)
    mu_1 <- mu_0 + c(X %*% alpha)
    PATE <- mean(X[s == 0,] %*% alpha)
  } else if (scenario == "exchange"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(s/(1 + exp(-X %*% lambda)) + (1 - s)/(1 + exp(-X %*% delta)))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(s*(X %*% beta0) + (1 - s)*(X %*% beta1))
    mu_1 <- mu_0 + c(X %*% alpha)
    PATE <- mean(X[s == 0,] %*% alpha)
  } else if (scenario == "ps-mis") {
    f_X <- c(1/(1 + exp(-U %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(s/(1 + exp(-U %*% lambda)) + (1 - s)/(1 + exp(-U %*% delta)))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(s*(X %*% beta0) + (1 - s)*(X %*% beta1))
    mu_1 <- mu_0 + c(X %*% alpha)
    PATE <- mean(X[s == 0,] %*% alpha)
  } else if (scenario == "out-mis"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(s/(1 + exp(-X %*% lambda)) + (1 - s)/(1 + exp(-X %*% delta)))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(s*(U %*% beta0) + (1 - s)*(U %*% beta1))
    mu_1 <- mu_0 + c(U %*% alpha)
    PATE <- mean(U[s == 0,] %*% alpha)
  } else if (scenario == "instrument"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(s/(1 + exp(-X %*% lambda)) + (1 - s)/(1 + exp(-X %*% delta)))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(U %*% beta0)
    mu_1 <- mu_0 + c(U %*% alpha)
    PATE <- mean(U[s == 0,] %*% alpha)
  } else if (scenario == "prognostic"){
    f_X <- c(1/(1 + exp(-U %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(1/(1 + exp(-U %*% lambda)))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(s*(X %*% beta0) + (1 - s)*(X %*% beta1))
    mu_1 <- mu_0 + c(X %*% alpha)
    PATE <- mean(X[s == 0,] %*% alpha)
  } else if (scenario == "sample-overlap") {
    f_X <- c(1/(1 + exp(-X %*% (4*gamma))))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(1/(1 + exp(-X %*% lambda)))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(U %*% beta0)
    mu_1 <- mu_0 + c(U %*% alpha)
    PATE <- mean(U[s == 0,] %*% alpha)
  } else if (scenario == "treat-overlap"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(1/(1 + exp( -X %*% (4*lambda))))
    z <- rbinom(n, 1, e_X) # treatment
    mu_0 <- c(U %*% beta0)
    mu_1 <- mu_0 + c(U %*% alpha)
    PATE <- mean(U[s == 0,] %*% alpha)
  }
  
  # potential outcomes
  y_tmp0 <- rnorm(n, mu_0, sqrt(sig2))
  y_tmp1 <- rnorm(n, mu_1, sqrt(sig2))
  y_pot <- cbind(y_tmp0, y_tmp1)
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  # create simulation dataset
  sim <- list(y = y, z = z, s = s, X = X, U = U, PATE = PATE)
  
  return(sim)
  
}

# Fits the balancing weights using a variety of methods
simfit <- function(idx = 1, simDat, ...) {
  
  # print(idx)
  
  dat <- simDat[,idx]
  PATE <- dat$PATE
  
  Y <- dat$y
  Z <- dat$z
  S <- dat$s
  X <- dat$X
  
  # Calibration - Transport

  cal_t <- try( calibrate(S = S, Y = Y, Z = Z, X = X, fusion = FALSE), silent = TRUE )
  est_t <- ifelse(!inherits(cal_t, "try-error") & cal_t$estimate < 10 & cal_t$estimate > -15, cal_t$estimate, NA)
  var_t <- ifelse(!inherits(cal_t, "try-error") & cal_t$estimate < 10 & cal_t$estimate > -15, cal_t$variance, NA)
  
  # Calibration - Fusion
  
  cal_f <- try( calibrate(S = S, Y = Y, Z = Z, X = X, fusion = TRUE), silent = TRUE )
  est_f <- ifelse(!inherits(cal_f, "try-error") & cal_f$estimate < 10 & cal_f$estimate > -15, cal_f$estimate, NA)
  var_f <- ifelse(!inherits(cal_f, "try-error") & cal_f$estimate < 10 & cal_f$estimate > -15, cal_f$variance, NA)
  
  # TMLE - Transport

  tmle_t <- try( tmle(S = S, Y = Y, Z = Z, X = X), silent = TRUE )
  tmle_est_t <- ifelse(!inherits(tmle_t, "try-error") & tmle_t$estimate < 10 & tmle_t$estimate > -15, tmle_t$estimate, NA)
  tmle_var_t <- ifelse(!inherits(tmle_t, "try-error") & tmle_t$estimate < 10 & tmle_t$estimate > -15, tmle_t$variance, NA)
  
  # Augmented - Transport
  
  aug_t <- try( augment(S = S, Y = Y, Z = Z, X = X, fusion = FALSE), silent = TRUE )
  aug_est_t <- ifelse(!inherits(aug_t, "try-error") & aug_t$estimate < 10 & aug_t$estimate > -15, aug_t$estimate, NA)
  aug_var_t <- ifelse(!inherits(aug_t, "try-error") & aug_t$estimate < 10 & aug_t$estimate > -15, aug_t$variance, NA)
  
  # Augmented - Fusion
  
  aug_f <- try( augment(S = S, Y = Y, Z = Z, X = X, fusion = TRUE), silent = TRUE )
  aug_est_f <- ifelse(!inherits(aug_f, "try-error") & aug_f$estimate < 10 & aug_f$estimate > -15, aug_f$estimate, NA)
  aug_var_f <- ifelse(!inherits(aug_f, "try-error") & aug_f$estimate < 10 & aug_f$estimate > -15, aug_f$variance, NA)
  
  # Combine Results
  
  tau <- c(tmle_est_t,aug_est_t,est_t,aug_est_f,est_f)
  se <- c(sqrt(tmle_var_t),sqrt(aug_var_t),sqrt(var_t),sqrt(aug_var_f),sqrt(var_f))
  
  return(list(tau = tau, se = se, PATE = PATE))
  
}
