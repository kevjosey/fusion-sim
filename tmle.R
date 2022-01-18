
tmle <- function(S, Z, Y, X, fusion = FALSE) {

  star <- (Y - min(Y, na.rm = T))/(max(Y, na.rm = T) - min(Y, na.rm = T))
  
  idx <- which(apply(X, 2, var) != 0) # remove intercept
  dat_y <- data.frame(star = star, unname(X[,idx]))
  dat_s <- data.frame(S = S, unname(X[,idx]))
  dat_z <- data.frame(Z = Z, unname(X[,idx]))

  n <- nrow(X)
  n_0 <- sum(S == 0)
  n_1 <- sum(S == 1)
  
  samp <- predict(glm(S ~ ., data = dat_s, family = binomial(link = "logit")), 
                  newdata = dat_s, type = "response")
  treat <- predict(glm(Z ~ ., data = dat_z, family = binomial(link = "logit"), subset = S == 1), 
                   newdata = dat_z, type = "response")
  
  g0w <- ((1 - treat) * samp) / (1 - samp)
  g1w <- (treat * samp) / (1 - samp)
  h0w <- (1 - Z) * S / g0w
  h1w <- Z * S / g1w
  
  if (fusion) {
    fit_0 <- glm(star ~ ., data = dat_y, family = quasibinomial(link = "logit"), subset = Z == 0)
    fit_1 <- glm(star ~ ., data = dat_y, family = quasibinomial(link = "logit"), subset = Z == 1)
  } else {
    fit_0 <- glm(star ~ ., data = dat_y, family = quasibinomial(link = "logit"), subset = Z == 0 & S == 1)
    fit_1 <- glm(star ~ ., data = dat_y, family = quasibinomial(link = "logit"), subset = Z == 1 & S == 1)
  }
  
  # initial predictions
  mu_tmp <- cbind(predict(fit_0, newdata = dat_y, type = "response"),
                  predict(fit_1, newdata = dat_y, type = "response"))
  
  mu <- cbind(Z*mu_tmp[,2] + (1 - Z)*mu_tmp[,1], mu_tmp)
  
  epsilon <- coef(glm(star ~ -1 + offset(qlogis(mu[,1])) + h0w + h1w, 
                      family = quasibinomial(link = "logit"), subset = S == 1))
    
  # Update initial prediction.
  mu1 <- plogis(qlogis(mu) + cbind(epsilon[1] * h0w + epsilon[2] * h1w,
                                   epsilon[1]/g0w, epsilon[2]/g1w))
  
  # Scale to original
  mu_new <- mu*(max(Y, na.rm = TRUE) - min(Y, na.rm = TRUE)) + min(Y, na.rm = TRUE)
  mu1_new <- mu1*(max(Y, na.rm = TRUE) - min(Y, na.rm = TRUE)) + min(Y, na.rm = TRUE)
  
  # Estimate
  tmle_est <- mean(mu1_new[S==0,3] - mu1_new[S==0,2])

  # Get efficient influence curve
  eic <- c((h1w - h0w) * (Y - mu_new[,1]) + I(S == 0)*mu1_new[,3] - I(S == 0)*mu1_new[,2] - tmle_est)/mean(I(S == 0))
  
  return(list(estimate = tmle_est, variance = var(eic) / n))
  
}
