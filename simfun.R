########################################
## PURPOSE: Simulation Functions Code ##
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

cfit_target <- function(cmat, 
                 target,
                 base_weights,
                 coefs_init = NULL,
                 optim_ctrl = list(maxit = 500, reltol = 1e-10),
                 ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  fn <- match.fun("lagrange_target")
  
  # initialize coefs
  if (is.null(coefs_init))
    coefs_init <- rep(0, times = ncol(cmat)) 
  else if (length(coefs_init) != ncol(cmat))
    stop("length(coefs_init) != ncol(cmat)")
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  opt <- stats::optim(coefs_init, fn, method = "BFGS",
                      cmat = cmat,
                      base_weights = base_weights,
                      target = target,
                      control = optim_ctrl, ...)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  
  weights <- c( 1 + exp(-cmat %*% coefs) )
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit_target"
  return(out)
  
}

esteq_target <- function(X, Y, Z, S, weights, base_weights, target, tau) {
  
  eq1 <- S*(2*Z - 1)*base_weights*weights*X
  eq2 <- S*base_weights*(Z*weights*X - X)
  eq3 <- S*base_weights*X - (1 - S)*X
  eq4 <- (1 - S)*(X - target)
  eq5 <- S*base_weights*weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4, eq5) 
  return(eq)
  
}

lagrange_target <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum((base_weights)*(cmat %*% coefs - exp(-cmat %*% coefs)))
  out <- -temp + sum( target * coefs)
  return(out)
  
}

# Fits the balancing weights using a variety of methods
simfit <- function(idx = 1, simDat) {
  
  dat <- simDat[,idx]
  PATE <- dat$PATE
  
  # Calibration
  
  S <- dat$s
  Y <- dat$y
  Y1 <- dat$y[S == 1]
  Z <- dat$z
  Z1 <- dat$z[S == 1]
  X <- dat$X
  X1 <- dat$X[S == 1,]
  X0 <- dat$X[S == 0,]

  n_1 <- sum(S)
  n_0 <- sum(1 - S)
  m <- ncol(X)
  
  target <- n_1*colMeans(X0)
  fit_tmp <- cfit(cmat = X1, target = target, distance = "entropy")
  
  q <- fit_tmp$weights
  cmat <- cbind(as.matrix((2*Z1 - 1)*X1), as.matrix( Z1*X1 ))
  b <- c(rep(0, times = ncol(X1)), c(t(X1) %*% q))
  fit_sent <- try( cfit_target(cmat = cmat, target = b, base_weights = q), silent = TRUE )
  
  if (!inherits(fit_sent, "try-error")) {
    
    A <- cbind((2*Z - 1)*X, Z*X)
    p <- rep(1, n_1 + n_0)
    p[S == 1] <- fit_sent$weights
    q <- rep(1, n_1 + n_0)
    q[S == 1] <- fit_sent$base_weights
    coefs_1 <- fit_sent$coefs
    coefs_2 <- fit_tmp$coefs
    
    tau_sent <- sum(S*(2*Z - 1)*p*q*Y)/sum(S*p*q*Z)
    
    pweights <- as.vector( -exp(-A %*% coefs_1) )
    qweights <- as.vector( -exp(-X %*% coefs_2) )
    U <- matrix(0, ncol = 4*m, nrow = 4*m)
    v <- rep(0, times = 4*m + 1)
    meat <- matrix(0, ncol = 4*m + 1, nrow = 4*m + 1)
    
    for (i in 1:(n_1 + n_0)) {
      
      U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] + S[i]*q[i]*pweights[i] * (A[i,] %*% t(A[i,]))
      U[1:m,(2*m + 1):(3*m)] <- U[1:m,(2*m + 1):(3*m)] + S[i]*(2*Z[i] - 1)*qweights[i]*p[i] * (X[i,] %*% t(X[i,]))
      U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] + S[i]*Z[i]*qweights[i]*p[i] * (X[i,] %*% t(X[i,])) - qweights[i] * (X[i,] %*% t(X[i,]))
      U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] + S[i]*qweights[i]*X[i,] %*% t(X[i,])
      
      U[(2*m + 1):(3*m),(3*m + 1):(4*m)] <- U[(2*m + 1):(3*m),(3*m + 1):(4*m)] + S[i]*diag(-1, m, m)
      U[(3*m + 1):(4*m),(3*m + 1):(4*m)] <- U[(3*m + 1):(4*m),(3*m + 1):(4*m)] + (1 - S[i])*diag(-1, m, m)
      
      v[1:(2*m)] <- v[1:(2*m)] + S[i]*(2*Z[i] - 1)*q[i]*pweights[i] * (Y[i] - Z[i]*tau_sent) * A[i,]
      v[(2*m+1):(3*m)] <- v[(2*m + 1):(3*m)] + S[i]*(2*Z[i] - 1) * qweights[i]*p[i] * (Y[i] - Z[i]*tau_sent) * X[i,]
      v[4*m + 1] <- v[4*m + 1] - S[i]*q[i]*p[i]*Z[i]
      meat <- meat + tcrossprod(esteq_target(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], weights = p[i],
                                             base_weights = q[i], target = target/n_1, tau = tau_sent))
      
      
    }
    
    invbread <- matrix(0, nrow = 4*m + 1, ncol = 4*m + 1)
    invbread[1:(4*m),1:(4*m)] <- U
    invbread[4*m + 1,] <- v
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error"))
      bread <- try(MASS::ginv(invbread))
    
    if (inherits(bread, "try-error")) {
      
      var_sent[j] <- NA
      cp_sent[j] <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread) 
      var_sent <- sandwich[4*m + 1, 4*m + 1]
      cp_sent <- as.numeric(tau_sent - sqrt(var_sent)*1.96 <= PATE & tau_sent + sqrt(var_sent)*1.96 >= PATE)
      
    }
    
  } else {
    
    tau_sent[j] <- NA
    var_sent[j] <- NA
    cp_sent[j] <- NA
    
  }
  
  # calib 
  
  entfit <- calib(X = X1, Z = Z1, target = target, simple = FALSE)
  entest <- try( calest_pate(obj = entfit, X = X, Y = Y, Z = Z, S = S), silent = TRUE )
  
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
  
  # AIPW
  mod_1 <- lm(Y1 ~ ., data = cdat[Z1 == 1,])
  m_1 <- predict(mod_1, cdat)
  mod_0 <- lm(Y1 ~ ., data = cdat[Z1 == 0,])
  m_0 <- predict(mod_0, cdat)        
  AIPW_1 <- sum(Z1*wts*(Y1 - m_1))
  AIPW_0 <- sum((1-Z1)*wts*(Y1 - m_0))
  norm_1 <- (sum(Z1*wts))^-1
  norm_0 <- (sum((1-Z1)*wts))^-1
  
  # GLM Results
  
  design <- svydesign(ids = ~ 1, weights = ~ wts, data = data.frame(wts = wts, Y1 = Y1, Z1 = Z1))
  smod <- svyglm(Y1 ~ Z1, design = design, family = gaussian)
  glmrslt <- coef(smod)[2]
  
  # Calib Result
  
  if (inherits(entest, "try-error")){
    
    entrslt <- NA
    entvar <- NA
    entcp <- NA
    
  } else { 
    
    entrslt <- entest$tau
    entvar <- entest$variance
    entcp <- as.numeric(entrslt - sqrt(entvar)*1.96 <= PATE & entrslt + sqrt(entvar)*1.96 >= PATE)
    
  }
  
  # Outcome Results
  
  outrslt <- mean(g_1 - g_0)
  
  # AIPW Results
  
  arg_1 <- norm_1*AIPW_1 + mean(g_1)
  arg_0 <- norm_0*AIPW_0 + mean(g_0)
  aipwrslt <- arg_1 - arg_0
  
  # Combine Results
  
  tau <- c(glmrslt, outrslt, aipwrslt, entrslt, tau_sent)
  cp <- c(entcp, cp_sent)
  
  return(list(tau = tau, cp = cp, PATE = PATE))
  
}
