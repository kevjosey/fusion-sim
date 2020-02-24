library(sandwich)
library(survey)

source("D:/Github/target-sim/target.R")
source("D:/Github/target-sim/tmle.R")
source("D:/Github/target-sim/simfun.R")
source("D:/Github/cbal/R/cbalance.R")


iter <- 500
n <- 1000
sig2 <- 5
y_scen <- "a"
z_scen <- "a"
s_scen <- "a"

set.seed(06261992)

simDat <- replicate(iter, hte_data(n = n, sig2 = sig2, y_scen = y_scen, z_scen = z_scen, s_scen = s_scen))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat)

tau_tmp <- do.call(rbind, estList[1,])
cp_tmp <- do.call(rbind, estList[2,])
PATE_tmp <- do.call(c, estList[3,])
colnames(tau_tmp) <- c("GLM", "OUT", "TMLE", "CAL", "FUS")

tau <- colMeans(tau_tmp, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
cp <- colMeans(cp_tmp, na.rm = TRUE)
PATE <- mean(PATE_tmp)
