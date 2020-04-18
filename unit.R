library(sandwich)
library(survey)

source("D:/Github/fusion-sim/transport.R")
source("D:/Github/fusion-sim/tmle.R")
source("D:/Github/fusion-sim/simfun.R")

iter <- 1000
n <- 1000
# n_1 <- 1000
# n_0 <- 1000
sig2 <- 20
y_scen <- "a"
z_scen <- "b"
s_scen <- "b"

set.seed(06261987)

# simDat <- replicate(iter, gen_data(n_1 = n_1, n_0 = n_0, sig2 = sig2, y_scen = y_scen, z_scen = z_scen))
simDat <- replicate(iter, hte_data(n = n, sig2 = sig2, y_scen = y_scen, z_scen = z_scen, s_scen = s_scen))

mean(do.call(c, simDat[6,]))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat)

tau_tmp <- do.call(rbind, estList[1,])
se_tmp <- do.call(rbind, estList[2,])
PATE_tmp <- do.call(c, estList[3,])
colnames(tau_tmp) <- c("TMLE-T", "ACAL-T", "CAL-T", "ACAL-F", "CAL-F")

tau <- colMeans(tau_tmp, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
se <- colMeans(se_tmp, na.rm = TRUE)
PATE <- mean(PATE_tmp)

colMeans(tau_tmp - PATE, na.rm = T)

cp_tmp <- matrix(NA, nrow = iter, ncol = 2)
cp_tmp[,1] <- as.numeric(tau_tmp[,3] - se_tmp[,1]*1.96 <= PATE & tau_tmp[,4] + se_tmp[,1]*1.96 >= PATE)
cp_tmp[,2] <- as.numeric(tau_tmp[,5] - se_tmp[,2]*1.96 <= PATE & tau_tmp[,5] + se_tmp[,2]*1.96 >= PATE)
cp <- colMeans(cp_tmp, na.rm = TRUE)
