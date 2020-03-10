library(sandwich)
library(survey)

source("~/Github/combine-sim/transport.R")
source("~/Github/combine-sim/tmle.R")
source("~/Github/combine-sim/simfun.R")

iter <- 500
n <- 10000
sig2 <- 5
y_scen <- "b"
z_scen <- "a"
s_scen <- "a"

set.seed(09301987)

simDat <- replicate(iter, hte_data(n = n, sig2 = sig2, y_scen = y_scen, z_scen = z_scen, s_scen))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat)

tau_tmp <- do.call(rbind, estList[1,])
se_tmp <- do.call(rbind, estList[2,])
PATE_tmp <- do.call(c, estList[3,])
colnames(tau_tmp) <- c("OUT", "TMLE", "CAL-T", "ACAL", "CAL-F")

tau <- colMeans(tau_tmp, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
se <- colMeans(se_tmp, na.rm = TRUE)
PATE <- mean(PATE_tmp)

cp_tmp <- matrix(NA, nrow = iter, ncol = 2)
cp_tmp[,1] <- as.numeric(tau_tmp[,3] - se_tmp[,1]*1.96 <= PATE & tau_tmp[,4] + se_tmp[,1]*1.96 >= PATE)
cp_tmp[,2] <- as.numeric(tau_tmp[,5] - se_tmp[,2]*1.96 <= PATE & tau_tmp[,5] + se_tmp[,2]*1.96 >= PATE)
cp <- colMeans(cp_tmp, na.rm = TRUE)
