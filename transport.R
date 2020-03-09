cfit <- function(cmat, target,
                 distance = c("entropy", "binary", "shifted", "transport"),
                 base_weights = NULL,
                 coefs_init = NULL,
                 optim_ctrl = list(maxit = 500, reltol = 1e-10),
                 ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  if (!(distance %in% c("entropy", "binary", "shifted", "transport")))
    stop("distance must be either \"entropy\", \"binary\", \"shifted\", or \"transport\"")
  
  if (distance == "binary") {
    fn <- match.fun(lagrange_bent)
  } else if (distance == "shifted") {
    fn <- match.fun(lagrange_sent)
  } else if (distance == "transport") {
    fn <- match.fun(lagrange_transport)
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
  } else if (distance == "transport") {
    weights <- c( 1 + exp(-cmat %*% coefs) )
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

esteq_transport <- function(X, Y, Z, S, weights, base_weights, n_1, n_0, tau) {
  
  eq1 <- S*(2*Z - 1)*base_weights*weights*X
  eq2 <- S*base_weights*(Z*weights*X - X)
  eq3 <- S*base_weights*X - (n_1/n_0)*(1 - S)*X
  eq4 <- S*base_weights*weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}

esteq_fusion <- function(X, Y, Z, S, weights, base_weights, n_1, n_0, tau) {
  
  eq1 <- (2*Z - 1)*base_weights*weights*X
  eq2 <- base_weights*(Z*weights*X - X)
  eq3 <- S*base_weights*X - (n_1/n_0)*(1 - S)*X
  eq4 <- base_weights*weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}

lagrange_ent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}

lagrange_bent <- function(coefs, cmat, target, base_weights) {
  
  weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(cmat %*% coefs)) )
  temp <- sum(weights*log(weights/base_weights) + (1 - weights)*log((1 - weights)/(1 - base_weights)))
  out <- -temp - sum(weights * cmat %*% coefs) + sum(target * coefs)
  return(out)
  
}

lagrange_sent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(cmat %*% coefs - (base_weights - 1)*exp(-cmat %*% coefs))
  out <- -temp + sum(target * coefs)
  return(out)
  
}

lagrange_transport <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum((base_weights)*(cmat %*% coefs - exp(-cmat %*% coefs)))
  out <- -temp + sum( target * coefs)
  return(out)
  
}
