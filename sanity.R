
Sig <- matrix(0, nrow = 2, ncol = 2)
diag(Sig) <- 1
mu_1 <- 10
mu_0 <- 8

out_1 <- out_2 <- rep(0, 1000)

for (i in 1:1000){

eval <- eigen(Sig, symmetric = TRUE)
y_init <- matrix(stats::rnorm(1000*2, 0, 1), nrow = 1000, ncol = 2) # iid potential outcomes
y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
y_pot <- y_tmp + cbind(rep(mu_0, 1000), rep(mu_1, 1000)) # include causal effect
z <- rbinom(1000, 1, 0.2)

y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
out_1[i] <- sum((1/mean(z))*z*y - 1/(mean(1 - z))*(1 - z)*y)/1000
out_2[i] <- sum(y_pot[,2] - y_pot[,1])/1000
}

mean(out_1)
var(out_1)
mean(out_2)
var(out_2)
