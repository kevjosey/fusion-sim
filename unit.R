library(sandwich)
library(survey)
library(SuperLearner)

source("~/Github/fusion-sim/transport.R")
source("~/Github/fusion-sim/tmle.R")
source("~/Github/fusion-sim/aug.R")
source("~/Github/fusion-sim/simfun.R")

iter <- 1000
n <- 1000
sig2 <- 4
scenario <- "out-exchange"

set.seed(42)

simDat <- replicate(iter, gen_data(n = n, sig2 = sig2, scenario = scenario))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat)

tau_tmp <- do.call(rbind, estList[1,])
se_tmp <- do.call(rbind, estList[2,])
PATE <- mean(do.call(c, estList[3,]))
colnames(tau_tmp) <- c("TMLE-T", "AUG-T", "CAL-T", "TMLE-F", "AUG-F", "CAL-F")
tau_tmp[is.infinite(tau_tmp) | is.nan(tau_tmp)] <- NA
se_tmp[is.infinite(se_tmp) | is.nan(se_tmp)] <- NA

tau <- colMeans(tau_tmp - PATE, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
se <- colMeans(se_tmp, na.rm = TRUE)

cp_tmp <- matrix(NA, nrow = iter, ncol = 6)
cp_tmp[,1] <- as.numeric(tau_tmp[,1] - se_tmp[,1]*1.96 <= PATE & tau_tmp[,1] + se_tmp[,1]*1.96 >= PATE)
cp_tmp[,2] <- as.numeric(tau_tmp[,2] - se_tmp[,2]*1.96 <= PATE & tau_tmp[,2] + se_tmp[,2]*1.96 >= PATE)
cp_tmp[,3] <- as.numeric(tau_tmp[,3] - se_tmp[,3]*1.96 <= PATE & tau_tmp[,3] + se_tmp[,3]*1.96 >= PATE)
cp_tmp[,4] <- as.numeric(tau_tmp[,4] - se_tmp[,4]*1.96 <= PATE & tau_tmp[,4] + se_tmp[,4]*1.96 >= PATE)
cp_tmp[,5] <- as.numeric(tau_tmp[,5] - se_tmp[,5]*1.96 <= PATE & tau_tmp[,5] + se_tmp[,5]*1.96 >= PATE)
cp_tmp[,6] <- as.numeric(tau_tmp[,6] - se_tmp[,6]*1.96 <= PATE & tau_tmp[,6] + se_tmp[,6]*1.96 >= PATE)
cp <- colMeans(cp_tmp, na.rm = TRUE)
names(tau) <- names(cp) <- colnames(tau_tmp)

tau
cp
apply(tau_tmp, 2, sd, na.rm = T)
