#######################################################################
### PURPOSE: Functions for Transportability using Entropy Balancing ###
### BY: Kevin Josey                                                 ###
#######################################################################

lagrange <- function(coefs, A, b) {
  
  temp <- sum(exp(-A %*% coefs))
  out <- temp + sum(b * coefs)
  return(out)
  
}

esteq <- function(X, Y, Z, weights, target, tau) {
  
  eq1 <- Z*weights*X - target
  eq2 <- (1 - Z)*weights*X - target
  eq3 <- Z*weights*(Y - tau) - (1 - Z)*weights*Y
  
  eq <- c(eq1, eq2, eq3) 
  return(eq)
  
}

esteq_simple <- function(X, Y, Z, weights, target, tau) {
  
  eq1 <- weights*X - target
  eq2 <- Z*weights*(Y - tau) - (1 - Z)*weights*Y
  
  eq <- c(eq1, eq2) 
  return(eq)
  
}

bootit <- function(Z) {

  tidx <- which(Z == 1)
  idx1 <- sample(tidx, size = tmin, replace = TRUE)
  
  cidx <- which(Z == 0)
  idx2 <- sample(cidx, size = tmin, replace = TRUE)
  
  idx_tmp <- c(idx1, idx2)
  
  idx3 <- sample.int(n, size = n - length(idx_tmp), replace = TRUE)
  idx <- c(idx_tmp, idx3)
  
}

calib <- function(X, Z, target, optim_ctrl = list(maxit = 500, reltol = 1e-10), complex = TRUE, ...) {
  
  if (!is.matrix(X))
    stop("X must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  fn <- match.fun(lagrange)
  
  if (complex) {
  
    A <- cbind(Z*X, (1 - Z)*X)
    b <- c(target, target)
    
  } else {
    
    A <- cbind(X, Z)
    b <- c(target, 1/2)
    
  }
  
  # initialize coefs
  coefs_init <- rep(0, times = ncol(A))
  opt <- stats::optim(coefs_init, fn, method = "BFGS", A = A,
                      b = b, control = optim_ctrl, ...)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  
  weights <- c( exp(-A %*% coefs) )
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              X = X, Z = Z, target = target,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl,
              complex = complex)
  
  class(out) <- "calib"
  return(out)
  
}

calest <- function(obj, Y, method = c("sandwich", "bootstrap"),boot_iter = 1000, ...) {
  
  if (!inherits(obj, "calib"))
    stop("obj must be of class \"calib\"")
  
  if (!(method %in% c("sandwich", "bootstrap")))
    stop("method must be either \"sandwich\" or \"bootstrap\"")

  weights <- obj$weights
  coefs <- obj$coefs
  X <- obj$X
  Z <- obj$Z
  tm <- obj$target
  n <- length(Z)
  m <- ncol(X)
  
  tau <- sum(weights*(2*Z - 1)*Y)/sum(Z*weights)
  sandie_conv <- TRUE
  
  if (method == "sandwich") {
    
    if(obj$complex) {
    
      U <- matrix(0, ncol = 2*m, nrow = 2*m)
      v <- rep(0, times = 2*m + 1)
      meat <- matrix(0, ncol = 2*m + 1, nrow = 2*m + 1)
  
      for (i in 1:n) {
  
        U[1:m,1:m] <- U[1:m,1:m] - Z[i] * weights[i] * (X[i,] %*% t(X[i,]))
        U[(m + 1):(2*m),(m + 1):(2*m)] <- U[(m + 1):(2*m),(m + 1):(2*m)] - (1 - Z[i]) * weights[i] * (X[i,] %*% t(X[i,]))
        v[1:m] <- v[1:m] - Z[i] * weights[i] * (Y[i] - tau) * X[i,]
        v[(m + 1):(2*m)] <- v[(m + 1):(2*m)] + (1 - Z[i]) * weights[i] * Y[i] * X[i,]
        v[2*m + 1] <- v[2*m + 1] - weights[i]*Z[i]
        s <- esteq(X = X[i,], Y = Y[i], Z = Z[i], weights = weights[i], target = tm/n, tau = tau)
        meat <- meat + s %*% t(s)
  
      }
  
      invbread <- matrix(0, nrow = 2*m + 1, ncol = 2*m + 1)
      invbread[1:(2*m),1:(2*m)] <- U
      invbread[2*m + 1, ] <- v
  
      bread <- try(solve(invbread), silent = TRUE)
  
      if (inherits(bread, "try-error")) {
        
        sandie_conv <- FALSE
        warning("Sandwich estimator is singular. Trying Bootstrap Estimator.")
        
      } else {
        
        sandwich <- bread %*% meat %*% t(bread)
        variance <- sandwich[2*m + 1, 2*m + 1]
        
      }
      
    } else {
      
      U <- matrix(0, ncol = m, nrow = m)
      v <- rep(0, times = m + 1)
      meat <- matrix(0, ncol = m + 1, nrow = m + 1)
      
      for (i in 1:n) {
        
        U[1:m,1:m] <- U[1:m,1:m] - weights[i] * (X[i,] %*% t(X[i,]))
        v[1:m] <- v[1:m] - (2*Z[i] - 1) * weights[i] * (Y[i] - Z[i]*tau) * X[i,]
        v[m + 1] <- v[m + 1] - weights[i]*Z[i]
        s <- esteq_simple(X = X[i,], Y = Y[i], Z = Z[i], weights = weights[i], target = tm/n, tau = tau)
        meat <- meat + s %*% t(s)
        
      }
      
      invbread <- matrix(0, nrow = m + 1, ncol = m + 1)
      invbread[1:m,1:m] <- U
      invbread[m + 1, ] <- v
      
      bread <- try(solve(invbread), silent = TRUE)
      
      if (inherits(bread, "try-error")) {
        
        sandie_conv <- FALSE
        warning("Sandwich estimator is singular. Trying Bootstrap Estimator.")
        
      } else {
        
        sandwich <- bread %*% meat %*% t(bread)
        variance <- sandwich[m + 1, m + 1]
        
      }
      
    }

  } 
  
  if (method == "bootstrap" | !(sandie_conv) ) {

    boot_iter <- floor(boot_iter)

    if (boot_iter <= 0 | length(boot_iter) != 1)
      stop("boot_iter must be a positive integer")

    boot_idx <- replicate(boot_iter, bootit(Z))
    boot_est <- apply( boot_idx, 2, function(idx, obj, Y, Z) {

      fit <- calib(X = obj$X[idx,], Z = obj$Z[idx],
                   target = obj$target,
                   coefs_init = obj$coefs,
                   optim_ctrl = obj$optim_ctrl, 
                   complex = obj$complex)

      rweights <- fit$weights
      est <- sum((2*Z[idx] - 1)*rweights*Y[idx])/sum(Z[idx]*rweights)
      return(est)

    }, obj = obj, Y = Y, Z = Z)

    variance <- var(boot_est)

  }
  
  out <- list(tau = tau, variance = variance)
  
}