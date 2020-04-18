cfit <- function(cmat, target,
                 distance = c("entropy", "binary", "shifted", "fusion"),
                 base_weights = NULL,
                 coefs_init = NULL,
                 optim_ctrl = list(maxit = 500, reltol = 1e-10),
                 ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  if (!(distance %in% c("entropy", "binary", "shifted", "fusion")))
    stop("distance must be either \"entropy\", \"binary\", \"shifted\", or \"fusion\"")
  
  if (distance == "binary") {
    fn <- match.fun(lagrange_bent)
  } else if (distance == "shifted") {
    fn <- match.fun(lagrange_sent)
  } else if (distance == "fusion") {
    fn <- match.fun(lagrange_fusion)
  } else { # distance == "entropy"
    fn <- match.fun(lagrange_ent)
  }
  
  if (is.null(base_weights)) { # initialize base_weights
    
    if (distance == "binary") {
      base_weights <- rep(1/2, nrow(cmat))
    } else if (distance == "shifted") {
      base_weights <- rep(2, nrow(cmat))
    } else { # distance == "entropy"
      base_weights <- rep(1, nrow(cmat))
    }
    
  } else if (length(base_weights) != nrow(cmat)) { 
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = ncol(cmat)) 
  } else if (length(coefs_init) != ncol(cmat)) {
    stop("length(coefs_init) != ncol(cmat)")
  }
  
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
  
  if (distance == "binary") {
    weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(cmat %*% coefs)) )
  } else if (distance == "shifted") {
    weights <- c( 1 + (base_weights - 1)*exp(-cmat %*% coefs) )
  } else if (distance == "fusion") {
    weights <- base_weights*c( 1 + exp(-cmat %*% coefs) )
  } else { # distance == "entropy"
    weights <- c( base_weights*exp(-cmat %*% coefs) )
  }
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              distance = distance,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}

esteq_transport <- function(S, X, Y, Z, p, base_weights, theta, tau) {
  
  eq1 <- S*(Z*p*X - theta)
  eq2 <- S*((1 - Z)*p*X - theta)
  eq3 <- (1 - S)*(base_weights*X - theta)
  eq4 <- S*p*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}

esteq_fusion <- function(X, Y, Z, S, p, q, base_weights, theta, tau) {
  
  eq1 <- p*q*Z*X - theta
  eq2 <- p*q*(1 - Z)*X - theta
  eq3 <- S*(q*X - theta)
  eq4 <- (1 - S)*(base_weights*X - theta)
  eq5 <- p*q*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4, eq5) 
  return(eq)
  
}

lagrange_ent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}

lagrange_fusion <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum((base_weights)*(cmat %*% coefs - exp(-cmat %*% coefs)))
  out <- -temp + sum( target * coefs)
  return(out)
  
}

# Estimation with individual-level data. Y and Z are arbitrary for the target sample, 
# but must be included. Set them to 0.
transport_estimate <- function(obj, S, X, Y1, Z1, ...) {
  
  if (!inherits(obj, "cfit"))
    stop("obj must be of class \"cfit\"")
  
  A <- obj$cmat
  weights <- obj$weights
  base_weights <- obj$base_weights
  coefs <- obj$coefs
  n_0 <- sum(1 - S)
  n_1 <- sum(S)
  n <- n_1 + n_0
  m <- ncol(X)
  theta <- obj$target[1:m]/n_1
  
  Y <- rep(1, times = n)
  Y[S == 1] <- Y1
  Z <- rep(1, times = n)
  Z[S == 1] <- Z1
  
  if (is.null(base_weights))
    base_weights <- rep(1, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")
  
  tau <- sum(S*(weights*(2*Z - 1)*Y)/sum(S*Z*weights))
    
  U <- matrix(0, ncol = 3*m, nrow = 3*m)
  v <- rep(0, times = 3*m + 1)
  meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)
  
  for (i in 1:n) {
    
    U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] - weights[i] * A[i,] %*% t(A[i,])
    
    U[1:m, (2*m + 1):(3*m)] <- U[1:m, (2*m + 1):(3*m)] - diag(S[i], m, m)
    U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] - diag(S[i], m, m)
    U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] - diag((1 - S[i]), m, m)
    
    v[1:(2*m)] <- v[1:(2*m)] - weights[i] * (2*Z[i] - 1) * (Y[i] - Z[i]*tau) * A[i,]
    v[3*m + 1] <- v[3*m + 1] - S[i]*weights[i]*Z[i]
    
    meat <- meat +  tcrossprod(esteq_transport(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                               p = weights[i], base_weights = base_weights[i], 
                                               theta = theta, tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = 3*m + 1, ncol = 3*m + 1)
  invbread[1:(3*m),1:(3*m)] <- U
  invbread[3*m + 1, ] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error")) {
    
    sandwich <- NA
    variance <- NA
    
  } else {
    
    sandwich <- bread %*% meat %*% t(bread)
    variance <- sandwich[3*m + 1, 3*m + 1]
    
  }
    
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}

# but must be included. Set them to 0.
fusion_estimate <- function(obj, base_obj, S, X, Y, Z, ...) {
  
  if (!inherits(obj, "cfit"))
    stop("obj must be of class \"cfit\"")
  
  A <- obj$cmat
  R <- base_obj$cmat
  weights <- obj$weights
  base_weights <- base_obj$base_weights
  
  q <- base_obj$weights
  p <- weights/q
  coefs_1 <- obj$coefs
  coefs_2 <- base_obj$coefs
  
  n_0 <- sum(1 - S)
  n_1 <- sum(S)
  n <- n_1 + n_0
  m <- ncol(X)
  theta <- obj$target[1:m]/n
  
  if (is.null(base_weights))
    base_weights <- rep(1, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")

  tau <- sum(weights*(2*Z - 1)*Y)/sum(weights*Z)
  
  dp <- as.vector( -exp(-A %*% coefs_1) )
  dq <- as.vector( -exp(-R %*% coefs_2) )
  U <- matrix(0, ncol = 4*m, nrow = 4*m)
  v <- rep(0, times = 4*m + 1)
  meat <- matrix(0, ncol = 4*m + 1, nrow = 4*m + 1)
  
  for (i in 1:(n_1 + n_0)) {
    
    U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] + q[i]*dp[i] * (A[i,] %*% t(A[i,]))
    
    U[1:(2*m),(2*m + 1):(3*m)] <- U[1:(2*m),(2*m + 1):(3*m)] + 
      dq[i]*p[i]*(A[i,] %*% t(R[i,]))
    U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] + 
      dq[i] * (R[i,] %*% t(R[i,]))
    
    U[1:m,(3*m + 1):(4*m)] <- U[1:m,(3*m + 1):(4*m)] - diag(1, m, m)
    U[(m + 1):(2*m),(3*m + 1):(4*m)] <- U[(m + 1):(2*m),(3*m + 1):(4*m)] - diag(1, m, m)
    U[(2*m + 1):(3*m),(3*m + 1):(4*m)] <- U[(2*m + 1):(3*m),(3*m + 1):(4*m)] - diag(S[i], m, m)
    U[(3*m + 1):(4*m), (3*m + 1):(4*m)] <- U[(3*m + 1):(4*m), (3*m + 1):(4*m)] - diag((1 - S[i]), m, m)
    
    v[1:(2*m)] <- v[1:(2*m)] + q[i]*dp[i] * (2*Z[i] - 1)*(Y[i] - Z[i]*tau)*A[i,]
    v[(2*m+1):(3*m)] <- v[(2*m + 1):(3*m)] + dq[i]*p[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau)*R[i,]
    v[4*m + 1] <- v[4*m + 1] - q[i]*p[i]*Z[i]

    meat <- meat + tcrossprod(esteq_fusion(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                           p = p[i], q = q[i], base_weights = base_weights[i],
                                           theta = theta, tau = tau))
    
    
  }
  
  invbread <- matrix(0, nrow = 4*m + 1, ncol = 4*m + 1)
  invbread[1:(4*m),1:(4*m)] <- U
  invbread[4*m + 1,] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    variance <- NA
  
  else {
    
    sandwich <- bread %*% meat %*% t(bread)
    variance <- sandwich[4*m + 1, 4*m + 1]
    
  }
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}
