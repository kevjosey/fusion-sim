#########################################
## TITLE: poc.R                        ##
## PURPOSE: Chapter 3 proof of concept ##
#########################################

library(cbal)

hte_data <- function(n, sig2, rho, y_scen = c("a", "b"), z_scen = c("a", "b")){
  
  # error variance
  R <- matrix(rho, nrow = 2, ncol = 2)
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
  
  # effect coefficients
  beta <- c(210, 27.4, 13.7, 13.7, 13.7)
  gamma <- c(20, -13.7, 0, 0, 13.7)
  
  # propensity score
  if (z_scen == "b") {
    e_X <- 1/(1 + exp( -(-u1 + 0.5*u2 - 0.25*u3 - 0.1*u4 ) ) )
  } else { # z_scen == "a"
    e_X <- 1/(1 + exp( -(-x1 + 0.5*x2 - 0.25*x3 - 0.1*x4 ) ) )
  }
  
  z <- rbinom(n, 1, e_X)
  
  if (y_scen == "b") {
    X <- cbind(rep(1, times = n), u1, u2, u3, u4)
  } else { # y_scen == "b"
    X <- cbind(rep(1, times = n), x1, x2, x3, x4)
  }
  
  # outcome mean
  mu_0 <- X%*%beta
  mu_1 <- X%*%(beta + gamma)
  
  tau <- t(apply(X[z == 1,], 2, mean)) %*% gamma
  
  # potential outcomes
  eval <- eigen(Sig, symmetric = TRUE)
  y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
  y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
  y_pot <- y_tmp + cbind(mu_0, mu_1) # include causal effect
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  # create simulation dataset
  sim <- list(y = y, z = z, x1 = x1, x2 = x2, x3 = x3, x4 = x4, ps = e_X)
  
  return(sim)
  
}

lagrange_ext <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum((base_weights)*(cmat %*% coefs - exp(-cmat %*% coefs)))
  out <- -temp + sum(target * coefs)
  return(out)
  
}

esteq_ext <- function(X, Y, Z, weights, base_weights, target, tau) {
  
  eq1 <- (2*Z - 1)*base_weights*weights*X
  eq2 <- base_weights*(Z*weights*X - X)
  eq3 <- base_weights*X - target
  eq4 <- base_weights*(Z*weights*(Y - tau) - (1 - Z)*weights*Y)
  
  eq <- c(eq1, eq2, eq3, eq4)
  return(eq)
  
}

set.seed(07271989)

tau <- 20
sig2 <- 10
rho <- 0
n <- 10000
iter <- 100

# simulate array of data
simDat <- replicate(iter, hte_data(n = n, sig2 = sig2, rho = rho, y_scen = "a", z_scen = "a")) # explores differences in ATE vs. ATT.

# fit using HTE specification
tau_sent <- var_sent <- cp_sent <- vector(mode = "numeric", length = iter)

for (j in 1:iter) {
  
  dat <- as.data.frame(simDat[,j])
  Y <- dat$y
  Z <- dat$z
  X <- model.matrix(z ~ x1 + x2 + x3 + x4, data = dat)
  n <- nrow(X)
  m <- ncol(X)
  
  b <- n*c(1, 0.5, 0, 0, 0.5)
  fit_tmp <- cfit(cmat = X, target = b, distance = "entropy")
  
  q <- fit_tmp$weights
  cmat <- cbind(as.matrix((2*Z - 1)*X), as.matrix( Z*X ))
  target <- c(rep(0, times = m), c(t(X) %*% q))
  fit_sent <- cfit(cmat = cmat, target = target, distance = "shifted", base_weights = q, test_var = TRUE)
  
  if (fit_sent$converged) {
    
    A <- cbind((2*Z - 1)*X, Z*X)
    p <- fit_sent$weights
    q <- fit_sent$base_weights
    coefs_1 <- fit_sent$coefs
    coefs_2 <- fit_tmp$coefs
    
    tau_sent[j] <- sum((2*Z -1)*p*q*Y)/sum(p*q*Z)
    
    pweights <- as.vector( -exp(-A %*% coefs_1) )
    qweights <- as.vector( -exp(-X %*% coefs_2) )
    U <- matrix(0, ncol = 3*m, nrow = 3*m)
    v <- rep(0, times = 3*m + 1)
    meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)
    
    for (i in 1:n) {
      
      U[1:m,1:(2*m)] <- U[1:m,1:(2*m)] + (2*Z[i] - 1) * q[i]*pweights[i] * (X[i,] %*% t(A[i,]))
      U[1:m,(2*m + 1):(3*m)] <- U[1:m,(2*m + 1):(3*m)] + (2*Z[i] - 1) * qweights[i]*p[i] * (X[i,] %*% t(X[i,]))
      
      U[(m+1):(2*m),1:(2*m)] <- U[(m+1):(2*m),1:(2*m)] + Z[i] * q[i]*pweights[i] * (X[i,] %*% t(A[i,]))
      U[(m+1):(2*m),(2*m+1):(3*m)] <- U[(m+1):(2*m),(2*m+1):(3*m)] + Z[i] * qweights[i]*p[i] * (X[i,] %*% t(X[i,])) -
        qweights[i]*(X[i,] %*% t(X[i,]))
      
      U[(2*m+1):(3*m),(2*m+1):(3*m)] <- U[(2*m+1):(3*m),(2*m+1):(3*m)] + qweights[i] * (X[i,] %*% t(X[i,]))
      
      v[1:(2*m)] <- v[1:(2*m)] + (2*Z[i] - 1) * q[i]*pweights[i] * (Y[i] - Z[i]*tau) * A[i,]
      v[(2*m+1):(3*m)] <- v[(2*m+1):(3*m)] + (2*Z[i] - 1) * qweights[i]*p[i] * (Y[i] - Z[i]*tau) * X[i,]
      v[3*m + 1] <- v[3*m + 1] - q[i]*p[i]*Z[i]
      meat <- meat + tcrossprod(esteq_ext(X = X[i,], Y = Y[i], Z = Z[i], weights = p[i], 
                                            base_weights = q[i], target = b/n, tau = tau_sent[j]))
      
    }
    
    invbread <- matrix(0, nrow = 3*m + 1, ncol = 3*m + 1)
    invbread[1:(3*m),1:(3*m)] <- U
    invbread[3*m + 1,] <- v
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error"))
      bread <- try(MASS::ginv(invbread))
    
    if (inherits(bread, "try-error"))
      stop("inversion of \"bread\" matrix failed")
    
    sandwich <- bread %*% meat %*% t(bread)
    var_sent[j] <- sandwich[3*m + 1, 3*m + 1]
    cp_sent[j] <- as.numeric(tau_sent[j] - sqrt(var_sent[j])*1.96 <= tau & tau_sent[j] + sqrt(var_sent[j])*1.96 >= tau)
    
  } else {
    
    tau_sent[j] <- NA
    var_sent[j] <- NA
    cp_sent[j] <- NA
    
  }
  
}

mean(tau_sent)
mean(cp_sent)