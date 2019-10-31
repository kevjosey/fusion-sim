library(sandwich)
library(survey)

source("D:/Dropbox (ColoradoTeam)/Projects/Transportability/Code/calib.R")
source("D:/Dropbox (ColoradoTeam)/Projects/Transportability/Code/simfun.R")

iter <- 100
n_0 <- 1000
n_1 <- 1000
prob <- 0.3
sig2 <- 5
scen <- "baseline"

set.seed(06261992)

simDat <- replicate(iter, gen_data(n_1 = n_1, n_0 = n_0, prob = prob, sig2 = sig2, scenario = scen))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat)

tau_tmp <- do.call(rbind, estList[1,])
cp_tmp <- do.call(rbind, estList[2,])
colnames(tau_tmp) <- c("GLM", "OUT", "AIPW", "MOM", "ENT")
colnames(cp_tmp) <- c("ENTSATE", "ENTPATE", "OUTSATE", "OUTPATE")

tau <- apply(tau_tmp, 2, mean, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
cp <- apply(cp_tmp, 2, mean, na.rm = TRUE)
