aug <- function(S, Z, Y, X, fusion = FALSE) {
  
  n <- nrow(X)
  n_0 <- sum(S == 0)
  n_1 <- sum(S == 1)
  theta <- colMeans(X[S == 0,])
  
  R <- as.matrix(S*X)
  fit_base <- cfit(cmat = R, target = n_1*theta, distance = "entropy")
  q <- fit_base$weights
  
  idx <- which(apply(X, 2, var) != 0) # remove intercept
  dat_y <- data.frame(Y = Y, unname(X[,idx]))
  dat_z <- data.frame(Z = Z, unname(X[,idx]))
  
  treat <- predict(glm(Z ~ ., data = dat_z, 
                       family = binomial(link = "logit"), 
                       subset = S == 1), 
                   newdata = dat_z, type = "response") 
  
  ipw <- ifelse(Z == 1, 1/treat, 1/(1-treat))
  
  if (fusion) {
    
    fit_0 <- lm(Y ~ ., data = dat_y, subset = Z == 0)
    fit_1 <- lm(Y ~ ., data = dat_y, subset = Z == 1)
    
    treat <- predict(glm(Z ~ ., data = dat_z, family = binomial(link = "logit")), 
                     newdata = dat_z, type = "response")
    ipw <- ifelse(Z == 1, 1/treat, 1/(1-treat))
    h0w <- (1 - Z) * (q*ipw)
    h1w <- Z * (q*ipw)
    
  } else {
    
    treat <- predict(glm(Z ~ ., data = dat_z, family = binomial(link = "logit"), subset = S == 1), 
                     newdata = dat_z, type = "response")
    ipw <- ifelse(Z == 1, 1/treat, 1/(1-treat))
    fit_0 <- lm(Y ~ ., data = dat_y, subset = Z == 0 & S == 1)
    fit_1 <- lm(Y ~ ., data = dat_y, subset = Z == 1 & S == 1)
    h0w <- (1 - Z) * S * (q*ipw)
    h1w <- Z * S * (q*ipw)
    
  }
  
  mu_tmp <- cbind(predict(fit_0, newdata = dat_y, type = "response"),
                  predict(fit_1, newdata = dat_y, type = "response"))
  
  mu <- cbind(Z*mu_tmp[,2] + (1 - Z)*mu_tmp[,1], mu_tmp)
  
  if (fusion) {
    aug_est <- sum((h1w - h0w) * (Y - mu[,1]))/sum(Z*q*ipw) + mean(mu[S == 0,3] - mu[S == 0,2])
    eic <- c((h1w - h0w) * (Y - mu[,1]) + (I(S == 0)*mu[,3] - I(S == 0)*mu[,2])/mean(I(S == 0)) - aug_est)
  } else {
    aug_est <- sum((h1w - h0w) * (Y - mu[,1]))/sum(S*Z*q*ipw) + mean(mu[S == 0,3] - mu[S == 0,2])
    eic <- c((h1w - h0w) * (Y - mu[,1]) + I(S == 0)*mu[,3] - I(S == 0)*mu[,2] - aug_est)/mean(I(S == 0))
  }
  
  return(list(estimate = aug_est, variance = var(eic) / n))
  
}
