
tmle <- function(data, S, Z, Y, nsitemodel,  nzmodel, noutmodel) {
  
  n.dat <- nrow(data)
  
  # Calculate components of clever covariate
  glm_cps <- glm(formula = nsitemodel, data = data, family = "binomial")
  cps <- predict(glm_cps, type = "response")
  
  glm_cpz <- glm(formula = nzmodel, data = data, family = "binomial")
  cpz <- predict(glm_cpz, type = "response")
  
  # Calculate clever covariate.
  ps0 <- mean(I(S == 0))
  
  g0w <- ((1 - cpz) * cps) / (1 - cps)
  g1w <- (cpz * cps) / (1 - cps)
  
  h0w <- ((1 - Z) * I(S == 1)) / g0w
  h1w <- (Z * I(S == 1)) / g1w 
  
  ymodel <- glm(formula = noutmodel, family = "gaussian", data = data, subset = S == 1)
  
  data_new0 <- data_new1 <- as.data.frame(data)
  data_new0$Z <- 0
  data_new1$Z <- 1
  
  # Initial prediction.
  q <- cbind(predict(ymodel, type = "link", newdata = data),
             predict(ymodel, type = "link", newdata = data_new0),
             predict(ymodel, type = "link", newdata = data_new1))
  
  epsilon <- coef(glm(Y ~ -1 + offset(q[, 1]) + h0w + h1w, family = "gaussian", subset = S == 1))
  
  # Update initial prediction.
  q1 <- q + c((epsilon[1] * h0w + epsilon[2] * h1w), epsilon[1] / g0w , epsilon[2] / g1w )
  
  # Get efficient influence curve values for everyone
  tmleest <- mean(q1[, 3][S == 0]) - mean(q1[, 2][S == 0])
  eic <- (((Z * h1w / ps0) - ((1 - Z) * h0w / ps0)) * (Y - plogis(q[, 1]))) + 
    (I(S == 0) / ps0 * plogis(q1[, 3])) - (I(S == 0) / ps0 * plogis(q1[, 2])) - 
    (tmleest / ps0)
  
  return(list(estimate = tmleest, variance = var(eic) / n.dat))
  
}
