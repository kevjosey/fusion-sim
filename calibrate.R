calibrate <- function(S, Z, Y, X, fusion = FALSE) {
  
  Y1 <- Y[S == 1]
  Z1 <- Z[S == 1]
  X0 <- X[S == 0,]
  theta <- colMeans(X0)
  
  n_1 <- sum(S)
  n_0 <- sum(1 - S)
  n <- n_1 + n_0
  m <- ncol(X)
  
  if (fusion) {
  
    A <- cbind(as.matrix(S*Z*X), as.matrix(S*(1 - Z)*X), as.matrix((1 - S)*Z*X), as.matrix((1 - S)*(1 - Z)*X))
    b <- c(n_1*theta, n_1*theta, n_0*theta, n_0*theta)
    fit_f <- cfit(cmat = A, target = b)
    est_f <- if(fit_f$converged) {
      fusion_estimate(obj = fit_f, theta = theta, S = S, X = X, Z = Z, Y = Y) 
    } else { list(estimate = NA, variance = NA) }
    
    return(est_f)
    
  } else {
    
    A <- cbind(as.matrix(S*(2*Z - 1)*X), as.matrix(S*X))
    b <- c(rep(0,m), n_1*theta)
    fit_t <- cfit(cmat = A, target = b)
    est_t <- if(fit_t$converged) {
      transport_estimate(obj = fit_t, theta = theta, S = S, X = X, Z1 = Z1, Y1 = Y1)
    } else{ list(estimate = NA, variance = NA) }
    
    return(est_t)
    
  }
  
}

cfit <- function(cmat, target,
                 base_weights = NULL,
                 coefs_init = NULL,
                 optim_ctrl = list(maxit = 500, reltol = 1e-10),
                 ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  fn <- match.fun(lagrange_ent)
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1, nrow(cmat))
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
  weights <- c( base_weights*exp(-cmat %*% coefs) )
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}

esteq_transport <- function(S, X, Y, Z, p, base_weights, theta, tau) {
  
  eq1 <- S*(2*Z - 1)*p*X
  eq2 <- S*(p*X - theta)
  eq3 <- (1 - S)*(base_weights*X - theta)
  eq4 <- S*p*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}

esteq_fusion <- function(S, X, Y, Z, p, base_weights, theta, tau) {
  
  eq1 <- S*(Z*p*X - theta)
  eq2 <- S*((1 - Z)*p*X - theta)
  eq3 <- (1 - S)*(Z*p*X - theta)
  eq4 <- (1 - S)*((1 - Z)*p*X - theta)
  eq5 <- (1 - S)*(base_weights*X - theta)
  eq6 <- p*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4, eq5, eq6) 
  return(eq)
  
}

lagrange_ent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}

# Estimation with individual-level data. Y and Z are arbitrary for the target sample, 
# but must be included. Set them to 0.
transport_estimate <- function(obj, theta, S, X, Y1, Z1, ...) {
  
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
  
  Y <- rep(1, times = n)
  Y[S == 1] <- Y1
  Z <- rep(1, times = n)
  Z[S == 1] <- Z1
  
  if (is.null(base_weights))
    base_weights <- rep(1, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")
  
  tau <- sum(S*weights*(2*Z - 1)*Y)/sum(S*Z*weights)
    
  U <- matrix(0, ncol = 3*m, nrow = 3*m)
  v <- rep(0, times = 3*m + 1)
  meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)
  
  for (i in 1:n) {
    
    U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] - weights[i] * A[i,] %*% t(A[i,])
    U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] - diag(S[i], m, m)
    U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] - diag(1 - S[i], m, m)
    
    v[1:(2*m)] <- v[1:(2*m)] - weights[i]*S[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau)*A[i,]
    v[3*m + 1] <- v[3*m + 1] - weights[i]*S[i]*Z[i]
    
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
fusion_estimate <- function(obj, theta, S, X, Y, Z, ...) {
  
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
  
  if (is.null(base_weights))
    base_weights <- rep(1, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")

  tau <- sum(weights*(2*Z - 1)*Y)/sum(Z*weights)
  
  U <- matrix(0, ncol = 5*m, nrow = 5*m)
  v <- rep(0, times = 5*m + 1)
  meat <- matrix(0, ncol = 5*m + 1, nrow = 5*m + 1)
  
  for (i in 1:n) {
    
    U[1:(4*m),1:(4*m)] <- U[1:(4*m),1:(4*m)] - weights[i] * A[i,] %*% t(A[i,])
    U[1:m,(4*m + 1):(5*m)] <- U[1:m,(4*m + 1):(5*m)] - diag(S[i], m, m)
    U[(m + 1):(2*m),(4*m + 1):(5*m)] <- U[(m + 1):(2*m),(4*m + 1):(5*m)] - diag(S[i], m, m)
    U[(2*m + 1):(3*m),(4*m + 1):(5*m)] <- U[(2*m + 1):(3*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    U[(3*m + 1):(4*m),(4*m + 1):(5*m)] <- U[(3*m + 1):(4*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    U[(4*m + 1):(5*m),(4*m + 1):(5*m)] <- U[(4*m + 1):(5*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    
    v[1:(4*m)] <- v[1:(4*m)] - weights[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau)*A[i,]
    v[5*m + 1] <- v[5*m + 1] - weights[i]*Z[i]
    
    meat <- meat +  tcrossprod(esteq_fusion(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                            p = weights[i], base_weights = base_weights[i], 
                                            theta = theta, tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = 5*m + 1, ncol = 5*m + 1)
  invbread[1:(5*m),1:(5*m)] <- U
  invbread[5*m + 1, ] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error")) {
    
    sandwich <- NA
    variance <- NA
    
  } else {
    
    sandwich <- bread %*% meat %*% t(bread)
    variance <- sandwich[5*m + 1, 5*m + 1]
    
  }
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}
