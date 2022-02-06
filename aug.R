aug <- function(S, Z, Y, X, fusion = FALSE) {
  
  idx <- which(apply(X, 2, var) != 0) # remove intercept
  dat_z <- data.frame(Z = Z, unname(X[,idx]))
  dat_y <- data.frame(Y = Y, unname(X[,idx]))
  
  n <- nrow(X)
  n_0 <- sum(S == 0)
  n_1 <- sum(S == 1)
  theta <- colMeans(X[S == 0,])
  
  # sampling scores
  R <- as.matrix(S*X)
  fit_base <- cfit(cmat = R, target = n_0*theta)
  q <- fit_base$weights

  # propensity scores
  treat <- predict(glm(Z ~ ., data = dat_z, family = binomial(link = "logit"), 
                       subset = S == 1), newdata = dat_z, type = "response")
  
  # outcome models
  if (fusion) {
    fit_0 <- lm(Y ~ ., data = dat_y, subset = Z == 0 & S == 0)
    fit_1 <- lm(Y ~ ., data = dat_y, subset = Z == 1 & S == 0)
  } else {
    fit_0 <- lm(Y ~ ., data = dat_y, subset = Z == 0 & S == 1)
    fit_1 <- lm(Y ~ ., data = dat_y, subset = Z == 1 & S == 1)
  }
  
  g0w <- (1 - treat) / q 
  g1w <- treat / q 
  h0w <- (1 - Z) * S / g0w
  h1w <- Z * S / g1w
  
  mu_tmp <- cbind(predict(fit_0, newdata = dat_y, type = "response"),
                  predict(fit_1, newdata = dat_y, type = "response"))
  
  mu <- cbind(Z*mu_tmp[,2] + (1 - Z)*mu_tmp[,1], mu_tmp)
  
  aug_est <- sum((h1w - h0w) * (Y - mu[,1]))/n_0 + mean(mu[S == 0,3] - mu[S == 0,2])
  eic <- c((h1w - h0w) * (Y - mu[,1]) + I(S == 0)*(mu[,3] - mu[,2] - aug_est))/mean(I(S == 0))
  
  return(list(estimate = aug_est, variance = var(eic) / n))
  
}
