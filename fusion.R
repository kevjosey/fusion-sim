cfit_fusion <- function(cmat, 
                           target,
                           base_weights,
                           coefs_init = NULL,
                           optim_ctrl = list(maxit = 500, reltol = 1e-10),
                           ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  fn <- match.fun("lagrange_fusion")
  
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
  
  tm <- target*nrow(cmat)
  
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
  
  class(out) <- "cfit_fusion"
  return(out)
  
}

esteq_fusion <- function(X, Y, Z, S, weights, base_weights, target, tau) {
  
  eq1 <- (2*Z - 1)*base_weights*weights*X
  eq2 <- base_weights*(Z*weights*X - X)
  eq3 <- S*base_weights*X - (1 - S)*X
  eq5 <- S*base_weights*weights*(Z*(Y - tau) - (1 - Z)*Y) + (1 - S)*weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4, eq5) 
  return(eq)
  
}

lagrange_fusion <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum((base_weights)*(cmat %*% coefs - exp(-cmat %*% coefs)))
  out <- -temp + sum( target * coefs)
  return(out)
  
}